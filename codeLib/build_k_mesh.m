%This function builds a mesh in fourier space from a mesh defined in real
%space, using MATLAB's convention for ordering of frequencies.

%nx - number of mesh columns
%ny - number of mesh rows
%dx - mesh element width in x-direction

function [kx, ky, dkx, dky] = build_k_mesh(nx,ny,dx)

dkx = 2*pi/(nx*dx);
dky = 2*pi/(ny*dx);

% kxv = [0:floor(nx/2), -ceil(nx/2)+1:-1]*dkx;
% kyv = [0:floor(ny/2), -ceil(ny/2)+1:-1]*dky;
% [kx, ky] = meshgrid(kxv,kyv);

kxv = [0:ceil(nx/2)-1, -floor(nx/2):-1]*dkx;
kyv = [0:ceil(ny/2)-1, -floor(ny/2):-1]*dky;
[kx, ky] = meshgrid(kxv,kyv);


end