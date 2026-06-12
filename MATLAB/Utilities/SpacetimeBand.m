function SpacetimeBand(P, Z, y_center, width)
%SpacetimeBand Creates spacetime diagrams by averaging over a y-band
%
%   SpacetimeBand(P, Z, y_center, width)
%
%   Inputs:
%       P           - 3D array (nx, ny, nt) of P values
%       Z           - 3D array (nx, ny, nt) of Z values
%       y_center    - Center of the band in y-direction [0,1]
%       width       - Width of the band in y-direction [0,1]

arguments
    P (:,:,:) double {mustBeReal}
    Z (:,:,:) double {mustBeReal}
    y_center (1,1) double {mustBeInRange(y_center, 0, 1)}
    width (1,1) double {mustBePositive, mustBeLessThanOrEqual(width, 1)}
end

% Check that P and Z have the same dimensions
if ~isequal(size(P), size(Z))
    error("P and Z must have the same dimensions");
end

[nx, ny, nt] = size(P);

% Physical coordinates
x = linspace(0,1,nx);
y = linspace(0,1,ny);

% Band limits in physical space
y_min = max(0, y_center - width/2);
y_max = min(1, y_center + width/2);

% Convert to indices
y0 = find(y >= y_min, 1, 'first');
y1 = find(y <= y_max, 1, 'last');

% Average over the selected band
spacetimeP = squeeze(mean(P(:, y0:y1, :), 2));
spacetimeZ = squeeze(mean(Z(:, y0:y1, :), 2));

figure();

ax1 = subplot(1,2,1);
imagesc(x, 1:nt, spacetimeP', 'Parent', ax1);
axis(ax1, 'xy');
colormap(ax1, parula);
colorbar(ax1);

xlabel(ax1, 'x');
ylabel(ax1, 'Time step');
title(ax1, sprintf('P (%.3f ≤ y ≤ %.3f)', y(y0), y(y1)));

ax2 = subplot(1,2,2);
imagesc(x, 1:nt, spacetimeZ', 'Parent', ax2);
axis(ax2, 'xy');
colormap(ax2, parula);
colorbar(ax2);

xlabel(ax2, 'x');
ylabel(ax2, 'Time step');
title(ax2, sprintf('Z (%.3f ≤ y ≤ %.3f)', y(y0), y(y1)));

end