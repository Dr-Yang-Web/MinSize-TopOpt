function H = filter_2d(rmin, nelx, nely)
% This educational code was written by Bo Yang, Shuntian Jia, Liang Gao, and Hao Li
% in June 2026 at Huazhong University of Science and Technology (Wuhan, China).
% Constructs a linear density filter for minimum length-scale control.
% Implements a circular convolution kernel with radius 'rmin'.

% Preallocation for sparse matrix construction
maxNZ = nelx*nely*(2*ceil(rmin)+1)^2;
iH = ones(maxNZ,1); jH = ones(maxNZ,1); sH = zeros(maxNZ,1);
k = 0;

% Build sparse filter matrix H
% Each row defines neighborhood weights for element e1
R = ceil(rmin);
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely + j1;
        for i2 = max(i1-R,1):min(i1+R,nelx)
            for j2 = max(j1-R,1):min(j1+R,nely)
                e2 = (i2-1)*nely + j2;
                k = k + 1;
                iH(k) = e1; jH(k) = e2;
                % Linear weight decay from center to radius rmin
                sH(k) = max(0, rmin - sqrt((i1-i2)^2 + (j1-j2)^2)) / rmin;
            end
        end
    end
end

% Assemble sparse matrix and normalize rows to sum to unity
H = sparse(iH(1:k), jH(1:k), sH(1:k));
Hs = sum(H, 2);
H = H ./ (Hs + (Hs == 0)); % Avoid division by zero
end