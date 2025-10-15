function Hb = Hb_estimator(RBC_per_kgbw)
% Hb_estimator Estimate haemoglobin (g/dL) from RBC per kg-bw (translated from Julia)
% Formula used in Julia:
%   nu_e = 8.2e-14; VL = 0.560; w = 7.0;
%   Hb = (100/3)*(nu_e / VL) * w * x;

nu_e = 8.2e-14;
VL = 0.560;
w = 7.0; % weight used in Julia example; if you want per-individual, pass weight separately
Hb = (100/3) * (nu_e / VL) * w .* RBC_per_kgbw;
end
