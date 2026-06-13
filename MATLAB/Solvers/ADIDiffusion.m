function [t, P, Z] = ADIDiffusion(f, g, P0, Z0, Dp, Dz, t_end, n)
% ADIDiffusion ADI method for reaction-diffusion problems.
%
%   [t, P, Z] = ADIDiffusion(f, g, P0, Z0, Dp, Dz, t_end, n)
%   Solves the problem using the ADI (Alternating Direction Implicit) method
%   for diffusion terms and explicit Euler for reaction terms.
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

% Build the 1D coefficient matrices for ADI scheme
% Using 1D Laplacian with Neumann BCs (2nd order)

% 1D Laplacian matrix for x-direction
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

% Identity matrices
Ix = speye(Nx);
Iy = speye(Ny);

% ADI matrices for P (x-direction implicit, y-direction explicit)
% and then x-direction explicit, y-direction implicit

% Matrices for x-direction implicit step
APx = Ix - (dt*Dp/2) * Lx;
% Matrices for y-direction implicit step
APy = Iy - (dt*Dp/2) * Ly;

% Matrices for Z
AZx = Ix - (dt*Dz/2) * Lx;
AZy = Iy - (dt*Dz/2) * Ly;

% Factorize matrices for each direction (LU decomposition)
% For P
[LPx, UPx] = lu(APx);
[LPy, UPy] = lu(APy);
% For Z
[LZx, UZx] = lu(AZx);
[LZy, UZy] = lu(AZy);

for i = 1:n-1
    Pk = P(:,:,i);
    Zk = Z(:,:,i);

    % Compute reaction terms at current time (explicit)
    R_P = f(t(i), Pk, Zk);
    R_Z = g(t(i), Pk, Zk);

    % ADI Step 1: x-direction implicit, y-direction explicit
    % Compute right-hand side for step 1 (y-direction explicit part)
    RHS1_P = zeros(Nx, Ny);
    RHS1_Z = zeros(Nx, Ny);

    for j = 1:Ny
        % For each y-line, apply the y-direction explicit operator
        tempP = Pk(:,j);
        tempZ = Zk(:,j);

        % Apply Ly to the column vector
        Ly_tempP = Ly * tempP;
        Ly_tempZ = Ly * tempZ;

        RHS1_P(:,j) = tempP + (dt*Dp/2) * Ly_tempP + dt * R_P(:,j);
        RHS1_Z(:,j) = tempZ + (dt*Dz/2) * Ly_tempZ + dt * R_Z(:,j);
    end

    % Solve x-direction implicit systems for each y-line
    P_star = zeros(Nx, Ny);
    Z_star = zeros(Nx, Ny);

    for j = 1:Ny
        P_star(:,j) = UPx \ (LPx \ RHS1_P(:,j));
        Z_star(:,j) = UZx \ (LZx \ RHS1_Z(:,j));
    end

    % ADI Step 2: y-direction implicit, x-direction explicit
    % Compute right-hand side for step 2 (x-direction explicit part)
    RHS2_P = zeros(Nx, Ny);
    RHS2_Z = zeros(Nx, Ny);

    for i_idx = 1:Nx
        % For each x-line, apply the x-direction explicit operator
        tempP = P_star(i_idx, :)';
        tempZ = Z_star(i_idx, :)';

        % Apply Lx to the row vector (converted to column)
        Lx_tempP = Lx * tempP;
        Lx_tempZ = Lx * tempZ;

        RHS2_P(i_idx, :) = (tempP + (dt*Dp/2) * Lx_tempP)';
        RHS2_Z(i_idx, :) = (tempZ + (dt*Dz/2) * Lx_tempZ)';
    end

    % Solve y-direction implicit systems for each x-line
    P_next = zeros(Nx, Ny);
    Z_next = zeros(Nx, Ny);

    for i_idx = 1:Nx
        P_next(i_idx, :) = (UPy \ (LPy \ RHS2_P(i_idx, :)'))';
        Z_next(i_idx, :) = (UZy \ (LZy \ RHS2_Z(i_idx, :)'))';
    end

    P(:,:,i+1) = P_next;
    Z(:,:,i+1) = Z_next;
end

end