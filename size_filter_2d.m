function H = size_filter_2d(rphys, nelx, nely,h)
% This educational code was written by Bo Yang, Shuntian Jia, Liang Gao, and Hao Li
% in June 2026 at Huazhong University of Science and Technology (Wuhan, China).
% Constructs a filter weights matrix for 2D minimum size control.
%   Adapted from the top88 code developed by Ole Sigmund’s team,
%   rphys: physical control radius; h: element size.
%   nelx, nely: #elements in x/y directions;.
%   ns: #elements inside ideal test geometry (circle intersection).
%   H(i,j)=1/ns if element j is within test circle of element i.

rtest = rphys/h; R = ceil(rtest); nele = nelx*nely; k = 0; % Initialization
maxNZ = nele*(2*R+1)^2; iH = ones(maxNZ,1); jH = iH; sH = zeros(maxNZ,1);

% Build sparse filter matrix H (neighbor weights for element e1)
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely + j1;
        for i2 = max(i1-R,1):min(i1+R,nelx)
            for j2 = max(j1-R,1):min(j1+R,nely)
                if rtest - norm([i1-i2,j1-j2]) >= 0
                    k = k+1; iH(k)=e1;
                    jH(k)=(i2-1)*nely+j2; sH(k)=1;
                end
            end
        end
    end
end

% Compute normalization factor (ns) based on forward cone volume
ns = 0;
for i = -R:R
    for j = -R:R
        if norm([i,j])<=rtest && j>=0 && norm([i,j-rtest])<=rtest, ns=ns+1; end
    end
end

H = sparse(iH,jH,sH) / ns; % Normalize filter weights