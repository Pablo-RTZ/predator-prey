function DiffusionPlot(t,P,Z)
%DiffusionPlot Animates the diffusion of P and Z over time
%
%   DiffusionPlot(t, P, Z)
%
%   Inputs:
%       t   - 1D array of time points (length n)
%       P   - 3D array (Nx, Ny, n) of P values over time
%       Z   - 3D array (Nx, Ny, n) of Z values over time

arguments
    t (:,1) double {mustBeReal, mustBeNonnegative}
    P (:,:,:) double {mustBeReal}
    Z (:,:,:) double {mustBeReal}
end

[Nx, Ny, n] = size(P);

% Fix color limits so colors don't rescale every frame
Pmin = min(P(:));
Pmax = max(P(:));

Zmin = min(Z(:));
Zmax = max(Z(:));

figure('Color','w');
x = linspace(0,1,Nx); y = linspace(0,1,Ny);

% Prey field
ax1 = subplot(1,2,1);
hP = imagesc(x, y, P(:,:,1)', 'Parent', ax1);
axis(ax1, 'image');
set(ax1, 'YDir', 'normal');
clim(ax1, [Pmin Pmax]);
colorbar(ax1);
xlabel(ax1, 'x');
ylabel(ax1, 'y');
titleP = title(ax1, sprintf('P(x,y,t), t = %.2f', t(1)));

% Predator field
ax2 = subplot(1,2,2);
hZ = imagesc(x, y, Z(:,:,1)', 'Parent', ax2);
axis(ax2, 'image');
set(ax2, 'YDir', 'normal');
clim(ax2, [Zmin Zmax]);
colorbar(ax2);
xlabel(ax2, 'x');
ylabel(ax2, 'y');
titleZ = title(ax2, sprintf('Z(x,y,t), t = %.2f', t(1)));

for k = 1:5:n
    % Update image data only
    hP.CData = P(:,:,k)';
    hZ.CData = Z(:,:,k)';

    % Update titles only
    titleP.String = sprintf('P(x,y,t), t = %.2f', t(k));
    titleZ.String = sprintf('Z(x,y,t), t = %.2f', t(k));

    drawnow
end