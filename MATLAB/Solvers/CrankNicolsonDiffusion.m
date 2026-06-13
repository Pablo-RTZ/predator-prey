function [t, P, Z] = CrankNicolsonDiffusion(f, g, P0, Z0, Dp, Dz, t_end, n)
% CrankNicolsonDiffusion Crank-Nicolson method for reaction-diffusion problems.
%
%   [t, P, Z] = CrankNicolsonDiffusion(f, g, P0, Z0, Dp, Dz, t_end, n)
%   Solves the problem using the Crank-Nicolson method for diffusion terms
%   and explicit Euler for reaction terms.
%
%   Inputs:
%       f       - Function handle for the prey reaction
%       g       - Function handle for the predator reaction
%       P0      - Initial condition for P
%       Z0      - Initial condition for Z
%       Dp      - Diffusion coefficient for P
%       Dz      - Diffusion coefficient for Z
%       t_end   - Time interval end
%       n       - Number of timesteps

arguments
    f (1,1) function_handle
    g (1,1) function_handle
    P0 (:,:) double
    Z0 (:,:) double
    Dp (1,1) double
    Dz (1,1) double
    t_end (1,1) double
    n (1,1) double {mustBeInteger,mustBePositive}
end

if ~isequal(size(P0), size(Z0))
    error('Initial conditions must have the same size.');
end

P(:,:,1) = P0;
Z(:,:,1) = Z0;
[Nx, Ny] = size(P0);

dx = 1/(Nx-1); dy = 1/(Ny-1);
x = linspace(0,1,Nx); y = linspace(0,1,Ny);

t = linspace(0, t_end, n);
dt = t_end/(n-1);

% Build the coefficient matrices for the implicit diffusion terms
% Using 5-point stencil for Laplacian with Neumann BCs
Ix = speye(Nx); Iy = speye(Ny);

% 1D Laplacian matrix with Neumann BCs (2nd order)
Lx = spalloc(Nx, Nx, 3*Nx);
Lx(1,1) = -2/dx^2; Lx(1,2) = 2/dx^2;
for i = 2:Nx-1
    Lx(i,i-1) = 1/dx^2;
    Lx(i,i) = -2/dx^2;
    Lx(i,i+1) = 1/dx^2;
end
Lx(Nx,Nx-1) = 2/dx^2; Lx(Nx,Nx) = -2/dx^2;

% 1D Laplacian matrix for y-direction
Ly = spalloc(Ny, Ny, 3*Ny);
Ly(1,1) = -2/dy^2; Ly(1,2) = 2/dy^2;
for i = 2:Ny-1
    Ly(i,i-1) = 1/dy^2;
    Ly(i,i) = -2/dy^2;
    Ly(i,i+1) = 1/dy^2;
end
Ly(Ny,Ny-1) = 2/dy^2; Ly(Ny,Ny) = -2/dy^2;

% 2D Laplacian using Kronecker product: L = Iy ⊗ Lx + Ly ⊗ Ix
L = kron(Iy, Lx) + kron(Ly, Ix);

% Crank-Nicolson matrices for diffusion terms
Mp = speye(Nx*Ny) - (dt*Dp/2) * L;
Mz = speye(Nx*Ny) - (dt*Dz/2) * L;

% Factorize matrices for efficiency (LU decomposition)
[Lp, Up] = lu(Mp);
[Lz, Uz] = lu(Mz);

for i = 1:n-1
    Pk = P(:,:,i);
    Zk = Z(:,:,i);
    
    % Compute reaction terms at current time (explicit)
    R_P = f(t(i), Pk, Zk);
    R_Z = g(t(i), Pk, Zk);
    
    % Compute diffusion terms at current time (explicit part of CN)
    P_vec = Pk(:);
    Z_vec = Zk(:);
    
    % Laplacian at current time
    LP_vec = L * P_vec;
    LZ_vec = L * Z_vec;
    
    % Right-hand side for Crank-Nicolson
    RHS_P = P_vec + (dt*Dp/2) * LP_vec + dt * R_P(:);
    RHS_Z = Z_vec + (dt*Dz/2) * LZ_vec + dt * R_Z(:);
    
    % Solve the linear systems
    P_next_vec = Up \ (Lp \ RHS_P);
    Z_next_vec = Uz \ (Lz \ RHS_Z);
    
    % Reshape back to 2D
    P(:,:,i+1) = reshape(P_next_vec, [Nx, Ny]);
    Z(:,:,i+1) = reshape(Z_next_vec, [Nx, Ny]);
end

end