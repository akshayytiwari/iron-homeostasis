function dydt = iron_odes(t, y, params)
% iron_odes  ODEs for paediatric iron homeostasis (translated from model_and_functions_file.jl)
% y order:
% 1 HMO, 2 FOS, 3 Inulin, 4 VitC, 5 LfFe, 6 LumenFe3, 7 LumenFe2, 8 LIP,
% 9 Ft, 10 FtFe, 11 B, 12 P, 13 RBC, 14 TFe, 15 H, 16 MFe

% unpack states
HMO     = y(1); FOS   = y(2); Inulin = y(3); VitC   = y(4);
LfFe    = y(5); LFe3  = y(6); LFe2   = y(7); LIP    = y(8);
Ft      = y(9); FtFe  = y(10); B      = y(11); P      = y(12);
RBC     = y(13); TFe  = y(14); H      = y(15); MFe    = y(16);

% for safety, if params provided as vector convert or expect struct fields
% We expect params to be a struct with fields named as in the Julia p[1]..p[56] mapping
% Example required fields (names used below) include:
% D_H, kappa_H, K_H, phi_H, D_F, kappa_F, K_F, phi_F, D_I, kappa_I, K_I, phi_I,
% D_V, phi_V, D_L, k_L, phi_L, D_Fe3, c_max, K_pH, K_Vit, phi_Fe3, D_Fe2,
% kappa_Fe, K_Fe, k_enter, beta, phi_Fe2, k_1, k_2, k_utilize, k_fpn, theta,
% mu_B, eta_B0, a_H, a_F, a_I, a_BP, mu_P, eta_P0, a_Fe, p_max, epsilon, zeta, delta_R,
% f_p, gamma, delta_Tfe, h, K_T, K_p, delta_H, delta_M, omega, nu
% For the actual field names check the params struct made by import_params_from_excel or adjust names below.

% To keep code readable, give local names from params.
% (Use try/catch if some fields missing)
fld = @(n) params.(n);

D_H = fld('D_H'); kappa_H = fld('kappa_H'); K_H = fld('K_H'); phi_H = fld('phi_H');
D_F = fld('D_F'); kappa_F = fld('kappa_F'); K_F = fld('K_F'); phi_F = fld('phi_F');
D_I = fld('D_I'); kappa_I = fld('kappa_I'); K_I = fld('K_I'); phi_I = fld('phi_I');
D_V = fld('D_V'); phi_V = fld('phi_V');
D_L = fld('D_L'); k_L = fld('k_L'); phi_L = fld('phi_L');
D_Fe3 = fld('D_Fe3'); c_max = fld('c_max'); K_pH = fld('K_pH'); K_Vit = fld('K_Vit'); phi_Fe3 = fld('phi_Fe3');
D_Fe2 = fld('D_Fe2'); kappa_Fe = fld('kappa_Fe'); K_Fe = fld('K_Fe'); k_enter = fld('k_enter'); beta = fld('beta'); phi_Fe2 = fld('phi_Fe2');
k_1 = fld('k_1'); k_2 = fld('k_2'); k_utilize = fld('k_utilize'); k_fpn = fld('k_fpn'); theta = fld('theta');
mu_B = fld('mu_B'); eta_B0 = fld('eta_B0'); a_H = fld('a_H'); a_F = fld('a_F'); a_I = fld('a_I'); a_BP = fld('a_BP');
mu_P = fld('mu_P'); eta_P0 = fld('eta_P0'); a_Fe = fld('a_Fe');
p_max = fld('p_max'); epsilon = fld('epsilon'); zeta = fld('zeta'); delta_R = fld('delta_R');
f_p = fld('f_p'); gamma = fld('gamma'); delta_Tfe = fld('delta_Tfe');
h = fld('h'); K_T = fld('K_T'); K_p = fld('K_p'); delta_H = fld('delta_H');
delta_M = fld('delta_M'); omega = fld('omega'); nu = fld('nu');

% Derived constants (same algebra as Julia)
a_HB = 1/kappa_H;
a_HP = -6/(7*kappa_H);
aF = 1.5/kappa_F;
aI = 1/kappa_I;
aFe = 2/kappa_Fe;

eta_B = eta_B0 * (1 + a_HB*kappa_H*HMO/(K_H + HMO)) * (1 + aF*kappa_F*FOS/(K_F + FOS)) * (1 + aI*kappa_I*Inulin/(K_I + Inulin));
eta_P = eta_P0 * (1 + a_HP*kappa_H*HMO/(K_H + HMO)) * (1 + aFe*kappa_Fe*LFe2/(K_Fe + LFe2));
pRBC = p_max * (epsilon / (epsilon + RBC)) * (zeta / (zeta + TFe));
pH = omega - nu*log10(max(B,1));  % clip B to avoid log10(0)

% ODEs
dHMO = D_H - kappa_H*HMO*(B + P)/(K_H + HMO) - phi_H*HMO;
dFOS = D_F - kappa_F*FOS*B/(K_F + FOS) - phi_F*FOS;
dInulin = D_I - kappa_I*Inulin*B/(K_I + Inulin) - phi_I*Inulin;
dVitC = D_V - phi_V*VitC;
dLfFe = D_L - k_L*LfFe - phi_L*LfFe;

% Lumen Fe3 -> Fe2 conversion enhanced by VitC and depends on pH (same algebra as Julia)
conv_term = c_max * LFe3 * (K_pH/(K_pH + pH)) * (1 + (VitC/(K_Vit + VitC)));
dLFe3 = D_Fe3 - conv_term - phi_Fe3*LFe3;

dLFe2 = D_Fe2 + conv_term - kappa_Fe*P.*LFe2./(K_Fe + LFe2) - k_enter*LFe2 * beta/(beta + FtFe) - phi_Fe2*LFe2;

dLIP = k_enter*LFe2 * beta/(beta + FtFe) + k_L*LfFe + k_1*FtFe - k_2*Ft.*LIP - k_utilize*LIP - k_fpn*LIP*theta/(theta + H);
dFt = k_1*FtFe - k_2*LIP.*Ft;
dFtFe = k_2*LIP.*Ft - k_1*FtFe;

dB = mu_B*B .* (1 - B./eta_B) - a_BP*B.*P;
dP = mu_P*P .* (1 - P./eta_P) - a_BP*B.*P;

dRBC = pRBC - delta_R*RBC;

dTFe = k_fpn*LIP*theta/(theta + H) - pRBC*f_p*TFe + delta_R*RBC*f_p*TFe*gamma/(gamma + H) - delta_Tfe*TFe;

dH = h*(TFe/(K_T + TFe)) .* (1 - (pRBC)/(K_p + pRBC)) - delta_H*H;

dMFe = delta_R*RBC*f_p*TFe*(1 - (gamma)/(gamma + H)) - delta_M*MFe;

dydt = [dHMO; dFOS; dInulin; dVitC; dLfFe; dLFe3; dLFe2; dLIP; dFt; dFtFe; dB; dP; dRBC; dTFe; dH; dMFe];
end
