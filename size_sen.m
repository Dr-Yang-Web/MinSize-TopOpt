function [Ls,Lv,dLs,dLv] = size_sen(volfrac,xPhys,Hs,Hv,H,dxPhys,model_domain,undesign,singular)
% This matlab program was written by Bo Yang, Liang Gao, Hao Li, Ying Zhou
% in June 2026 at Huazhong University of Science and Technology (Wuhan, China).
% Computes minimum size constraint (Ls/Lv) and their design sensitivities (dLs/dLv)
% H: filter weights matrix; dxPhys: derivative of projection function in three-field SIMP.
% model_domain: Element position index of the model (no extension domain included)
% undesign/singular: Undesign-domain/singular-region element position index


kappa = 200; lamda = 0.95; % Projection steepness and threshold
ne = length(model_domain); % Number of elements in the geometric model

%% ----------------- SOLID-PHASE SIZE CONSTRAINT --------------------------
% Filter inverted physical variables and project
S = Hs*xPhys(:);
Theta = 0.5 + 0.5*tanh(kappa*(S-lamda));

% Enforce singular/undesign regions to satisfy size constraint
Theta(undesign) = 1; Theta(singular) = 1;

% Calculate solid-phase constraint function Ls (based on additive volume)
Gammas = (xPhys(:).*Theta).^2;
Ls = 1 - sum(Gammas(model_domain))/ne/volfrac;

%Sensitivity analysis to design variables: chain rule for dLs/dx
dLs0 = -1/ne/volfrac;
dTheta = (kappa/2) * sech(kappa * (S - lamda)).^2;
dGammas = 2*(xPhys(:).*Theta).*(Theta + Hs*(xPhys(:).*dTheta));
dLs = H*(dLs0.*dGammas.*dxPhys(:));

%% ----------------- VOID-PHASE/FILLET-RADII SIZE CONSTRAINT --------------
% Filter inverted physical variables and project
Sv = Hv*(1-xPhys(:));
Thetav = 0.5 + 0.5*tanh(kappa*(Sv-lamda));

% Calculate void-phase constraint function Lv
Gammav = ((1-xPhys(:)).*Thetav).^2;
Lv = 1 - sum(Gammav(model_domain))/ne/(1-volfrac);

% Sensitivity analysis to design variables: chain rule for dLv/dx
dLv0 = -1/ne/(1-volfrac);
dThetav = (kappa/2) * sech(kappa * (Sv - lamda)).^2;
dGammav = 2*((1-xPhys(:)).*Thetav).*(-Thetav + Hv*((1-xPhys(:)).*dThetav));
dLv = H*(dLv0.*dGammav.*dxPhys(:)); % Chain rule with H and dxPhys
end