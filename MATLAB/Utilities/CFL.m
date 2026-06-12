function [dt, n] = CFL(Dp,Dz,Nx,Ny,t_end)
%CFL Calculates the minimum timestep (and number of points) to assure
% stability for diffusion using explicit Euler and central differences
%
%   [dt, n] = CFL(Dp,Dz,Nx,Ny,t_end)
%
%   Inputs:
%       Dp  - Diffusion coefficient for P
%       Dz  - Diffusion coefficient for Z
%       Nx  - Number of nodes in x
%       Ny  - Number of nodes in y
%       t_end   - Time interval end

arguments
Dp (1,1) double {mustBePositive}
Dz (1,1) double {mustBePositive}
Nx (1,1) double {mustBeInteger,mustBePositive}
Ny (1,1) double {mustBeInteger,mustBePositive}
t_end (1,1) double {mustBePositive} = []
end

dx = 1/(Nx-1);
dy = 1/(Ny-1);

dt = 1/(2*max(Dp,Dz)*(1/dx^2+1/dy^2));

if ~isempty(t_end)
    n = ceil(t_end/dt);
else
    n = [];
end

end