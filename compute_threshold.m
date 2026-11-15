function t = compute_threshold(xg, sub)
% This matlab program was written by Bo Yang, Liang Gao, Hao Li, Ying Zhou
% in April 2026 at Huazhong University of Science and Technology (Wuhan, China).
% Computes isosurface threshold matching the target volume fraction.
% Uses sub-pixel interpolation and root-finding to ensure accuracy.

if nargin < 2, sub = 20; end  % Default subsampling factor

% Convert element-wise field to nodal field
xn = E2N(xg); [ny, nx] = size(xn);

% Create high-resolution grid for interpolation
x = linspace(0, nx-1, (nx-1)*sub) + 0.5/sub;
y = linspace(0, ny-1, (ny-1)*sub) + 0.5/sub;
[X, Y] = meshgrid(x, y);

% Interpolate density field to fine grid
V = interp2(1:nx, 1:ny, xn, X+1, Y+1);

% Define area fraction function and target volume
area_fun = @(t) mean(V(:) > t);
target = mean(xg(:));

% Solve for threshold where area fraction matches target
if target >= 1, t = 0;
elseif target <= 0, t = 1;
else
    t = fzero(@(t) area_fun(t) - target, 0.5);
end
end