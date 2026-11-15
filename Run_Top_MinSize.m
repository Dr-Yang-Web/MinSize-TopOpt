%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This educational code was written by Bo Yang, Shuntian Jia, Liang Gao, and Hao Li
% in June 2026 at Huazhong University of Science and Technology (Wuhan, China).
% The code is mainly adapted from top88 matlab code (www.topopt.dtu.dk.)
% developed by Ole Sigmund team in 2011.
% The MMA optimization solver is from Krister Svanberg in 2007 or 2008.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all; warning('off','all'); tic

%% Problem definition and parameters
folder = 'MBB_top_data_rs=4_rv=2_nx=600'; % Data saving directory
rs = 4; rv = 4; rf = rs;          % Solid/void min radius, filter radius
flag = 2;                         % Constraint flag: 0:none, 1:solid-only, 2:bi-phase
Wext = max(rs,rv);                % Domain extension width
nx = 600; ny = 100; ne = nx*ny;   % Design domain size
nelx = nx + 2*Wext; nely = ny + 2*Wext; nele = nelx*nely; % FE mesh size
penal = 3; Vsolid = 0.4;          % SIMP penalty and volume fraction
baita_max = 16; yita = 0.5;       % Projection steepness and threshold
E0 = 1; Emin = 1e-9; nu = 0.3;    % Material properties

%% PREPARE FINITE ELEMENT ANALYSIS
% Element stiffness matrix (4-node quad)
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

% Node and element connectivity
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

% Boundary conditions (MMBB beam)
F = zeros(2*(nely+1)*(nelx+1),1);
F(2*(nely+1)*nelx/2+2*(Wext+1)) = -1;% Mid-point load
fixeddofs1 = 2*(nely+1)*(Wext+1)-2*Wext+[-1 0]; % Bottom-left supports
fixeddofs2 = 2*(nely+1)*(nelx+1-Wext)-2*Wext; % Bottom-right support
fixeddofs = [fixeddofs1 fixeddofs2];
U = zeros(2*(nely+1)*(nelx+1),1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);

%% DEFINE DOMAIN INDICES (Passive, Undesign, Model)
% Passive region (fixed solid)
x0 = zeros(nely,nelx);
x0([1:Wext nely-Wext+1:nely],:) = 1;
x0(:,[1:Wext nelx-Wext+1:nelx]) = 1;
passive = find(x0==1);

% Undesign and singular regions
x1 = zeros(nely,nelx);
x1(nely-Wext-rs+1:nely-Wext,[Wext+1:Wext+rs nelx-Wext-rs+1:nelx-Wext]) = 1;
undesign = find(x1==1); singular = undesign;

% Active model domain
xd = zeros(nely,nelx);
xd(Wext+1:nely-Wext,Wext+1:nelx-Wext)=1;
model_domain = find(xd==1);

%% PREPARE FILTERS
H = filter_2d(rf, nelx, nely);      % Density filter
Hs = size_filter_2d(rs, nelx, nely,1); % Solid-phase size filter
Hv = size_filter_2d(rv, nelx, nely,1); % Void-phase size filter

%% INITIALIZATION
x = repmat(Vsolid, nely, nelx); % Design variables
maxloop = 300;                  % Max iterations
xsave = zeros(maxloop, nele); xPhys_save = zeros(maxloop,nele); % History storage
comp = zeros(1, maxloop); vol = zeros(1, maxloop);
loop = 0;

%% MMA optimizer setup
m = 1; n = nele; xmin = zeros(n,1); xmax = ones(n,1);
xval = x(:); xold1 = x(:); xold2 = x(:); a0 = 1;
low = zeros(n,1); upp = ones(n,1);

%% START ITERATION
while loop < 300
    loop = loop + 1;

    %% Parameter continuation (Projection steepness & constraint relaxation)
    baita = min(1+max(ceil((loop-30)/15),0),baita_max);
    sigma1 = max(2^(-floor(loop/15)),1e-4);  sigma2 = sigma1;

    %% Update three-field variables
    x(passive) = 0; x(undesign) = 1;
    xf = H*x(:); xf(passive) = 0; xf(undesign) = 1;
    xp = (tanh(baita*yita)+tanh(baita*(xf(:)-yita)))...
        ./(tanh(baita*yita)+tanh(baita*(1-yita)));
    xp(passive) = 0; xp(undesign) = 1;
    xPhys = reshape(xp(:),nely,nelx);   
    Mnd = 100 * sum(4*xp.*(1-xp)) / nele; % Measure of non-discreteness

    %% FE-ANALYSIS
    sK = reshape(KE(:)*((1-xPhys(:)').^penal*Emin ...
        + xPhys(:)'.^penal*E0),64*nelx*nely,1);
    K = sparse(iK,jK,sK);     K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
   
    %% OBJECTIVE AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum(((1-xPhys).^penal*Emin + xPhys.^penal*E0).*ce));

    dc0 = -(-penal*Emin*(1-xPhys).^(penal-1) + penal*E0*xPhys.^(penal-1)).*ce;
    dv0 = ones(nely,nelx);
    dp = baita*(1-tanh(baita*(yita-xf(:))).^2)...
        ./(tanh(baita*yita)+tanh(baita*(1-yita)));
    dp(passive) = 0;  dp(undesign) = 0;
    dc = H*(dc0(:).*dp(:));  dv = H*(dv0(:).*dp(:));

    % Size constraint functions and sensitivities
    [gms, gmv, dgs, dgv] = size_sen(Vsolid, xPhys, Hs,...
        Hv, H, dp, model_domain, undesign, singular);
 
    %% MMA SUBPROBLEM
    f0val = c; df0dx = dc(:);
    vt = sum(xPhys(model_domain))/ne;
    fval = vt/Vsolid - 1;
    dfdx = dv(:)'/ne/Vsolid;

    % Add size constraints when discreteness is sufficient
    if Mnd <= 20
        if flag == 1 % Solid-phase constraint
            m = 2;
            fval = [vt/Vsolid - 1; gms - sigma1];
            dfdx = [dv(:)'/ne/Vsolid; dgs(:)'];
        elseif flag == 2 % Bi-phase constraint
            m = 3;
            fval = [vt/Vsolid - 1; gms - sigma1; gmv - sigma2];
            dfdx = [dv(:)'/ne/Vsolid; dgs(:)'; dgv(:)'];
        end
    end
        
    a = zeros(m,1); c_MMA = 1e3*ones(m,1); d = zeros(m,1);
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low,upp] = ...
        mmasub(m, n, loop, xval, xmin, xmax, xold1, xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);

    %% UPDATE DESIGN VARIABLES
    xnew = reshape(xmma, nely, nelx);
    xold2 = xold1(:); xold1 = x(:);
    xnew(passive) = 0; xnew(undesign) = 1; 
    x = xnew; xval = x(:);
    change = sum(abs(xnew(:)-xold1(:)))/sum(abs(xnew(:)));

    %% Save history data
    xsave(loop,:) = x(:)'; xPhys_save(loop,:) = xPhys(:)';
    vol(loop) = vt; comp(loop) = c;

    %% Data visualization

    % Printing
    fprintf(' It.:%4i Obj.:%5.2f Vol.:%4.3f ch.:%5.4f Mnd.:%3.2f\n', ...
        loop, c, vt, change, Mnd);

    % Drawing
    drawnow;
    subplot(2,1,1);
    xs = reshape(xPhys(model_domain),ny,nx);
    colormap(gray); imagesc(1-xs); clim([0 1]); axis equal off;
    title(['Iteration: ', num2str(loop)]);

    subplot(2,1,2);
    [AX, H1, H2] = plotyy(1:loop, comp(1:loop), 1:loop, vol(1:loop));
    set(H1, 'LineWidth', 1); set(H2, 'LineWidth', 1); S = 12;
    set(AX, 'xlim', [-20 300]); set(AX(1), 'xtick', -20:40:300);
    set(AX(2), 'ylim', [0.3 0.5]);  set(AX(2), 'ytick', 0.3:0.1:0.5);
    title('Iteration process', 'FontSize', S);
    ylabel(AX(1), 'Objective', 'FontSize', S);
    ylabel(AX(2), 'Volume fraction', 'FontSize', S);
    legend([H1 H2], {'Objective', 'Volume fraction'}, ...
        'Location', 'best','FontSize', S-1);

    %% Convergence check
    if loop > 15 && mean(abs(comp(loop-15:loop)...
            -comp(loop-15))/comp(loop-15)) < 1e-5
        break;
    end
end

%% POST-PROCESSING: BINARY STRUCTURE AND OBJECTIVE
xd = sort(xp(model_domain), 'descend');
va = xd(round(ne*Vsolid)); % Threshold
xend = xPhys >= va; % Binarization
sK = reshape(KE(:) * ((1-xend(:)')*Emin + xend(:)'*E0), 64*nelx*nely, 1);
K = sparse(iK, jK, sK); K = (K + K')/2;
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
ce = reshape(sum((U(edofMat)*KE) .* U(edofMat), 2), nely, nelx);
c1 = sum(sum(((1-xend)*Emin + xend*E0) .* ce));
fprintf('True structural compliance: %5.2f\n', c1);

%% Smooth geometry representation via level-set
xg = reshape(xf(model_domain),ny,nx);
xn1 = E2N(xg); th = compute_threshold(xg); % Node density & threshold
figure; contourf(flipud(xn1), [th th]); colormap([0.7,0.7,0.7]); axis equal off;

%% Save data
mkdir(folder); sf = {'NC','SC','BC'}; sf = sf{flag+1};
eval(['xsave_' sf '=xsave; xPhys_' sf '=xPhys_save;' ...
    ' xcomp_' sf '=comp; xvol_' sf '=vol;']);
save(fullfile(folder, ['data_' sf '.mat']), ['xsave_' sf],...
    ['xPhys_' sf], ['xcomp_' sf], ['xvol_' sf]);