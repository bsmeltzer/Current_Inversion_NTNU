 
%clear 
close all;
addpath(genpath('C:\Sites\CopterCurrents'));
addpath(genpath('C:\Sites\Current_Inversion_NTNU'));


%% Loop over latitude bins and calculate k-dependent Doppler shift velocities

dataPath = 'C:\Sites\Current_Inversion_NTNU\exampleData';%Path to cube directory
cubeFileName = 'exampleCube.mat';%file name of .mat file with image cube.

clear STAT_SA;
verbose = 1;%Whether or not to show the surface images

load(fullfile(dataPath,cubeFileName));
fileName = RES.filename;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the IMG_SEQ structure. - step 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear IMG_SEQ;


for i = 1:size(RES.MAT,1)
    if verbose
        figure(1);imagesc(squeeze(RES.MAT(i,:,:)));colorbar;axis image;colormap gray;%caxis([0,150]);
    end
    %xlim([350,450]);ylim([300,500]);
drawnow;
IMG_SEQ.IMG(:,:,i) = DJIFrame_to_IMG_SEQ_format(squeeze(RES.MAT(i,:,:)));
end

[gridY,gridX] = meshgrid(1:size(IMG_SEQ.IMG,2),1:size(IMG_SEQ.IMG,1));
gridX = (gridX-mean(gridX(:)))*RES.dx;
gridY = (gridY-mean(gridY(:)))*RES.dy;

IMG_SEQ.gridX = gridX;%X-Y interchanging due to CopterCurrents convention
IMG_SEQ.gridY = gridY;%X-Y interchanging due to CopterCurrents convention
IMG_SEQ.dx = RES.dy;
IMG_SEQ.dy = RES.dx;
IMG_SEQ.dt = RES.dt;

%The rest of the fields should be unnecessary but populated with dummy variables
%to avoid causing errors.
IMG_SEQ.altitude = 100;
IMG_SEQ.pitch = -90;
IMG_SEQ.roll = 0;
IMG_SEQ.heading = 0;
IMG_SEQ.mediainfo_string = '';


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the wavenumber-dependent Doppler shift velocities - step 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = 1;%Whether or not to show diagnostic figures of SNR for each wavenumber.

dk = 2*pi/(IMG_SEQ.dx*min(size(IMG_SEQ.IMG,1),size(IMG_SEQ.IMG,2)));% wavenumber resolution of spectrum in each spatial window (not strictly true if dx ~= dy, but value only needs to be approximate in practice)
kW = 2*dk;%Half width of wavenumbers bins [rad/m]

%Wavenumber values for Doppler shift extraction.
%wavenumbers = dk*10:kW:maxFrequency^2/9.81;
wavenumbers = 6*dk:kW:1.6;
%wavenumbers = 1.0;
%wavenumbers = 6*kW:4*kW:2.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define parameters specific to the k-dependent Doppler shift extraction and
%PEDM

frequencyLimits = [0.5,6.0];%frequency limits for masking the spectrum [min max], rad/sec

% Ux current limits [m/s]
Ux_limits = 1.0*[-1.0 1.0];

% Uy current limts [m/s]
Uy_limits = 1.0*[-1.0 1.0];

% Current step size [m/s]
U_res = 0.1;

% (optional): whether to include 2nd harmonic of the spectrum in the fit (false by default)
include2ndHarmonic = [];

% (optional): whether to do the fit in log space (false by default)
logFlag = [];

% (optional) omegaWidthFun: function handle as a function of wavenumber i.e.
%@(k) f(k)...., specifying frequency width of the weighting function in
%frequency-angle space (constant wavenumber). Width is half-width 1/e^2
%point of a Gaussian function.
%omegaWidthFun = @(k) 0.4 + 0.1*k;
omegaWidthFun = @(k) 0.4 + 0.0*k;


%The following OPTIONAL parameters involve post-processing of the Doppler shifts:
%SNR_filter: whether to use a signal-to-noise filter (false by default)
SNR_filter = 0;

%SNR_threshold: threshold signal-to-noise value for above filter (set to 2.0 by default)
SNR_threshold = sqrt(1);

%Peak_filter: whether to use a multiple peaks filter (false by default)
Peak_filter = 0;

%Peak_threshold: peak threshold of maximum value (0.5 by default)
Peak_threshold = 0.5;

%Outlier_filter: whether to use an outlier filter (quartile-based) (false by default)
Outlier_filter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define the fluid & physical properties: depth, gravitational acceleration,
%surface tension coefficient
Properties = struct('h',Inf,'g',9.81,'T',0.072 / 1000);
[Uym, Uxm] = meshgrid(min(Uy_limits):U_res:max(Uy_limits),min(Ux_limits):U_res:max(Ux_limits));

%Initialize Doppler shift variables
Ux = zeros(1,numel(wavenumbers))*NaN;
Uy = zeros(1,numel(wavenumbers))*NaN;
UxC = zeros(1,numel(wavenumbers))*NaN;
UyC = zeros(1,numel(wavenumbers))*NaN;
UxLS = zeros(1,numel(wavenumbers))*NaN;
UyLS = zeros(1,numel(wavenumbers))*NaN;

SNR_max = zeros(1,numel(wavenumbers))*NaN;

clear DSV;
clear DSVPC;
clear DSVLS;

for jj = 1:numel(wavenumbers)
%%
wavenumberLimits = wavenumbers(jj) + kW*[-1,1];%wavenumber limits for masking the spectrum [min max], rad/sec

STCFIT = struct('Name',{fileName});
STCFIT = generate_STCFIT_for_NSPP(STCFIT,wavenumbers(jj),include2ndHarmonic,logFlag,...
    omegaWidthFun,SNR_filter,SNR_threshold,Peak_filter,Peak_threshold,Outlier_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now extract the Doppler shifts 
% 
% tic;
%STCFIT = run_current_fit_depth_profile(IMG_SEQ,STCFIT,windowList);
% toc;

% get spectrum
%Read in the power spectrum.

Spectrum_crop = retrieve_power_spectrum(...
     IMG_SEQ.IMG,....
     IMG_SEQ.dx,...
     IMG_SEQ.dy,...
     IMG_SEQ.dt,wavenumberLimits,frequencyLimits);
 


%%
fit_param = STCFIT.NSPP_fit_param;
fit_param.kWidth = 4*kW;

fit_param.Ux_2D = Uxm;
fit_param.Uy_2D = Uym;
%Get Doppler shift velocities
%tic;
out_DS = get_doppler_shift_velocities_nsp(Spectrum_crop,fit_param,Properties,verbose);

%toc;
DSV(jj) = out_DS;

Ux(jj) = -out_DS.Ux_filt;
Uy(jj) = -out_DS.Uy_filt;
SNR_max(jj) = out_DS.SNR_max;


dtheta = 2*dk / wavenumbers(jj);
Scyl = cylinder_cross_section(Spectrum_crop,dtheta,4);
thetaVals = Scyl.thetaM(1,:);


figure(5);imagesc(Scyl.thetaM(1,:)/pi,Scyl.omegaM(:,1),Scyl.P_k);colorbar;shading flat;set(gca,'YDir','normal');
line(thetaVals/pi,sqrt(9.81*wavenumbers(jj)) + 0*thetaVals,'LineStyle',':','Color','r','LineWidth',1.5);
line(thetaVals/pi,sqrt(9.81*wavenumbers(jj))+wavenumbers(jj)*(cos(thetaVals)*Ux(jj) + sin(thetaVals)*Uy(jj)),'LineStyle','--','Color','r','LineWidth',1.5);
%line(thetaVals/pi,sqrt(9.81*wavenumbers(jj))-wavenumbers(jj)*(cos(thetaVals)*UxLS(jj) + sin(thetaVals)*UyLS(jj)),'LineStyle','-','Color','r','LineWidth',1.5);

% hold on;plot(thetaVals,sqrt(9.81*wavenumbers(jj))-wavenumbers(jj)*(cos(thetaVals)*Ux(jj) + sin(thetaVals)*Uy(jj)),'--w',');
% hold off;
xlabel('\theta/\pi');ylabel('Frequency [rad/s]');
title(sprintf('$k$: %.2f Rad/m \n$\\Delta\\omega/k$: %.2f m/s',wavenumbers(jj),mean(diff(Scyl.omegaM(:,1))) / wavenumbers(jj)),'FontWeight','normal','Interpreter','latex');
drawnow;
%%


figure(8);
%yyaxis left;
plot(...
    wavenumbers,Ux,'o',...
    wavenumbers,Uy,'xk',...
    'MarkerSize',6);drawnow;
xlabel('Wavenumber [Rad/m]');ylabel('[m/s]');
ylim(1*[-1,1]);
% yyaxis right;
% plot(wavenumbers,SNR_max,'--');
% ylabel('SNR');
hl = legend('$\tilde{c}_x$','$\tilde{c}_y$','Location','sw');
%hl = legend('$\tilde{c}_x$','$\tilde{c}_y$','$\tilde{c}_x$','$\tilde{c}_y$','SNR','Location','e');

set(hl,'Interpreter','latex');
set(gca,'FontSize',14);
drawnow;





end


STAT_SA.DSV = DSV;
STAT_SA.nCubes = 1;

%save(fullfile(dataPath,'latbin_DopplerShifts_20220109.mat'),'STAT_SA');

%end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform the PEDM - step 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run inversion methods
inds = 1:numel(wavenumbers);
wavenumbers = [STAT_SA.DSV.wavenumbers];

Ux = -[STAT_SA.DSV.Ux];
Uy = -[STAT_SA.DSV.Uy];
out = find_current_depth_profile(wavenumbers(inds),Ux(inds),Uy(inds),Inf);
STAT_SA.PEDM = out;

%%% Step 7: plot depth profiles for an individual window.
[f1,f2,p1,p2] = plot_depth_profile_results(out,3);
%plot_pedm_results(out,4);



%% Wave spectrum

clear S_twk;
    
wnS = (wavenumbers(end)+kW)/dk * 2 + 1;
kSpec = zeros(wnS);
jj_skip = max(round(2*kW / mean(diff(wavenumbers))),1);

kxVec = -(wnS-1)/2*dk:dk:(wnS-1)/2*dk;
kyVec = -(wnS-1)/2*dk:dk:(wnS-1)/2*dk;

wni = 0;

for jj = 1:jj_skip:numel(wavenumbers)
%%
wavenumberLimits = wavenumbers(jj) + kW*[-1,1];%wavenumber limits for masking the spectrum [min max], rad/sec

STCFIT = struct('Name',{fileName});
STCFIT = generate_STCFIT_for_NSPP(STCFIT,wavenumbers(jj),include2ndHarmonic,logFlag,...
    omegaWidthFun,SNR_filter,SNR_threshold,Peak_filter,Peak_threshold,Outlier_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now extract the Doppler shifts 
% 
% tic;
%STCFIT = run_current_fit_depth_profile(IMG_SEQ,STCFIT,windowList);
% toc;

% get spectrum
Spectrum = retrieve_power_spectrum(...
    IMG_SEQ.IMG,....
    IMG_SEQ.dx,...
    IMG_SEQ.dy,...
    IMG_SEQ.dt,wavenumberLimits,frequencyLimits);

%%
fit_param = STCFIT.NSPP_fit_param;
fit_param.kWidth = 4*kW;

fit_param.Ux_2D = Uxm;
fit_param.Uy_2D = Uym;
%Get Doppler shift velocities
%tic;
% out_DS = get_doppler_shift_velocities_nsp(Spectrum,fit_param,Properties,verbose);
% 
% if savePlots
% SNR_fileName = fullfile(SNR_plotDir,...
%     sprintf('%04iDeg_%04iRadpkm_SNR',round(latitude*1e4),round(wavenumbers(jj)*1e3)));
% print(SNR_fileName,'-dpng');
% end
% %toc;
% DSV(jj) = out_DS;

FP = Properties;FP.omegaWidth = fit_param.omegaWidthFun(wavenumbers(jj));FP.kWidth = fit_param.kWidth;
[SNR,P_k,G,inds] = nsp_doppler_shift_extraction(Spectrum,FP,wavenumbers(jj),STAT_SA.DSV(jj).Ux,STAT_SA.DSV(jj).Uy);

WNS_j = zeros(size(Spectrum.power_Spectrum));
WNS_j(inds) = Spectrum.power_Spectrum(inds).*G;

WNSS = Spectrum;
WNSS.power_Spectrum = WNS_j;
dtheta = 2*dk / wavenumbers(jj);
Scyl = cylinder_cross_section(Spectrum,dtheta,4);
Scyl2 = cylinder_cross_section(WNSS,dtheta,4);

Scyl.k = wavenumbers(jj);
wni = wni + 1;
S_twk(wni) = Scyl;

WNS_j = nansum(WNS_j,3);
WNS_j = padarray(WNS_j,(wnS-size(WNS_j))/2,'both');
%size(Spectrum.power_Spectrum)

kSpec = kSpec + WNS_j;

figure(5);
subplot(2,1,1);imagesc(Scyl.thetaM(1,:)/pi,Scyl.omegaM(:,1),log10(Scyl.P_k));colorbar;shading flat;set(gca,'YDir','normal');
xlabel('\theta/\pi');ylabel('Frequency [rad/s]');
title(sprintf('$k$: %.2f Rad/m \n$\\Delta\\omega/k$: %.2f m/s',wavenumbers(jj),mean(diff(Scyl.omegaM(:,1))) / wavenumbers(jj)),'FontWeight','normal','Interpreter','latex');
drawnow;

subplot(2,1,2);imagesc(Scyl2.thetaM(1,:)/pi,Scyl2.omegaM(:,1),Scyl2.P_k);colorbar;shading flat;set(gca,'YDir','normal');
xlabel('\theta/\pi');ylabel('Frequency [rad/s]');
title(sprintf('$k$: %.2f Rad/m \n$\\Delta\\omega/k$: %.2f m/s',wavenumbers(jj),mean(diff(Scyl.omegaM(:,1))) / wavenumbers(jj)),'FontWeight','normal','Interpreter','latex');
drawnow;

figure(6);imagesc(-kxVec,kyVec,10*log10(kSpec));axis image;colorbar;
%caxis([-40,-5]);
xlabel('$k_y$ [rad/m]','Interpreter','latex');
ylabel('$k_x$ [rad/m]','Interpreter','latex');
drawnow;




end



kSpec2 = kSpec;
kSpec2(kSpec2 == 0) = min(kSpec2(kSpec2>0));
kSpec2 = imrotate(kSpec2,270);

%%

%figure(7);contourf(kxVec,-kyVec,imgaussfilt(10*log10(kSpec2),2),[-5:-2.5:-40]);axis image;colorbar;
figure(7);contourf(kxVec,-kyVec,imgaussfilt(10*log10(kSpec2),2));axis image;colorbar;
% 
% kxa = kxVec./abs(kxVec).*log10(abs(kxVec));kxa(isnan(kxa)) = 0;
% kya = -kyVec./abs(kyVec).*log10(abs(kyVec));kya(isnan(kya)) = 0;
% figure(7);contourf(...
%     kxa,kya,...
%     imgaussfilt(10*log10(kSpec2),2));axis image;colorbar;


xlabel('$k_x$ [rad/m]','Interpreter','latex');
ylabel('$k_y$ [rad/m]','Interpreter','latex');
grid on;
drawnow;


%% Plot the k-omega spectrum for a specific angle.



for thetaDeg = 20:-20:-170;

theta = thetaDeg*pi/180;
dtheta = 0.1*pi;

Stheta = zeros(size(S_twk(1).omegaM,1),numel(S_twk));

omegaVals = S_twk(1).omegaM(:,1);
kVals = [S_twk(:).k];

for i = 1:numel(S_twk)
    
    mask_i = abs(S_twk(i).thetaM-theta)<dtheta;
    S_i = S_twk(i).P_k.*mask_i;
    
    Stheta(:,i) = sum(S_i,2,'omitnan');
    
end

Ux = -[STAT_SA.DSV(:).Ux];
Uy = -[STAT_SA.DSV(:).Uy];

g = 9.8;
omega0 = sqrt(g*wavenumbers);
omegaDSP = omega0 + wavenumbers.*(cos(theta)*Ux + sin(theta)*Uy);
figure(8);pcolor(kVals,omegaVals,10*log10(Stheta));colorbar;shading flat;
hold on;
plot(wavenumbers,omegaDSP,'--r',...
    2*wavenumbers,2*omegaDSP,'-.r');
hold off;
xlabel('Wavenumber [rad/m]');ylabel('Frequency [rad/s]');
title(sprintf('Angle: %.0fDeg',round(thetaDeg)));
drawnow;

end



%% Plot the velocity shifts.

T = 30;
L = 128;
kvect = wavenumbers;

Cg = 0.5*sqrt(g./wavenumbers);
delta_k = 2*pi / L;
delta_om = 2*pi / T;

delC_k = Cg * delta_k./kvect;
delC_om = delta_om./kvect;

figure(1);loglog(kvect,delC_k,kvect,delC_om,'--');
hl = legend('$\Delta c_{\Delta k}$','$\Delta c_{\Delta\omega}$','Location','ne');
set(hl,'Interpreter','latex');
set(gca,'FontSize',12);
xlabel('Wavenumber [rad/m]');ylabel('[m/s]');
ylim([1e-2,1e1]);
title(sprintf('T = %i s, L = %i m',T,L));













