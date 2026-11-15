function nodal_density = E2N(element_density)
% This educational code was written by Bo Yang, Shuntian Jia, Liang Gao, and Hao Li
% in June 2026 at Huazhong University of Science and Technology (Wuhan, China).
% Extends element-wise density field to nodal density (2D/3D compatible).
% Input:  element_density - density at element centroids (nx×ny or nx×ny×nz)
% Output: nodal_density   - density at nodes (nx+1×ny+1 or nx+1×ny+1×nz+1)

switch ndims(element_density)
    case 2, nodal_density = expand_2d(element_density);
    case 3, nodal_density = expand_3d(element_density);
    otherwise, error('Input must be a 2D or 3D matrix.');
end
end

%% 2D NODAL DENSITY INTERPOLATION
function nodal = expand_2d(element)
[nx, ny] = size(element);
nodal = zeros(nx+1, ny+1);
[I, J] = meshgrid(1:nx+1, 1:ny+1);

% Assign average density from adjacent elements to each node
for idx = 1:numel(I)
    i = I(idx); j = J(idx);
    neighbors = element(max(i-1,1):min(i,nx), max(j-1,1):min(j,ny));
    nodal(i,j) = mean(neighbors(:));
end
end

%% 3D NODAL DENSITY INTERPOLATION
function nodal = expand_3d(element)
[nx, ny, nz] = size(element);
nodal = zeros(nx+1, ny+1, nz+1);
[I, J, K] = ndgrid(1:nx+1, 1:ny+1, 1:nz+1);

% Assign average density from adjacent elements to each node
for idx = 1:numel(I)
    i = I(idx); j = J(idx); k = K(idx);
    neighbors = element(max(i-1,1):min(i,nx), ...
                       max(j-1,1):min(j,ny), ...
                       max(k-1,1):min(k,nz));
    nodal(i,j,k) = mean(neighbors(:));
end
end