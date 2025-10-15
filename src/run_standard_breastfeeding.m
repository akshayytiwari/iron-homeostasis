% RUN_STANDARD_BREASTFEEDING  Simulate the model with periodic breastfeeding doses,
% mimicking the Julia master file behaviour.
%
% Place the parameter Excel/CSV into the same folder and set 'param_file' accordingly.

clear; clc; close all;

param_file = 'parameters_2mo_09Nov2024_supplement_foods_tuned.xlsx';
params = import_params_from_excel(param_file);

% set simulation / feeding schedule (these match the Julia master file defaults)
avg_weight = 5.0; % kg (example for 2-month)
freq = 12;        % feeds per day
n_days = 1000;    % long run to reach steady state

% compute dose amounts per feed (if params supply daily totals use those instead)
% The Julia master uses:
% tot_HMO = 8.0 (g/day) => dose_HMO = tot_HMO*1000/(avg_weight*freq) mg/kgbw per feed
tot_HMO = 8.0;
tot_LfFe = 0.16e-3; % g/day in Julia master
dose_HMO = tot_HMO*1000/(avg_weight*freq);
dose_LfFe = tot_LfFe*1000/(avg_weight*freq);

% initial conditions (taken from Julia master init_conds)
y0 = [0;0;0;0;0;0;0;0;8e-3;0;1e4;1e4;0;0;0;0];

% integrate between dose times and add doses manually (simple event handling)
dose_times = 0 : 1/freq : n_days;
t_all = [];
y_all = [];

t0 = 0;
y = y0;
options = odeset('RelTol',1e-6,'AbsTol',1e-9);

for i = 2:length(dose_times)
    tspan = [dose_times(i-1), dose_times(i)];
    [tseg,yseg] = ode15s(@(t,y) iron_odes(t,y,params), tspan, y, options);
    t_all = [t_all; tseg(1:end-1)];
    y_all = [y_all; yseg(1:end-1,:)];
    % apply dose at t = dose_times(i)
    y = yseg(end,:)';
    y(1) = y(1) + dose_HMO;  % HMO
    y(5) = y(5) + dose_LfFe; % LfFe
end

% append last point
t_all = [t_all; dose_times(end)];
y_all = [y_all; y'];

% save steady state vector (last column)
ss = y';
writematrix(ss, 'ssvals_standard_breastfeeding_matlab.csv');
fprintf('Saved steady-state vector to ssvals_standard_breastfeeding_matlab.csv\n');

% compute Hb from RBC steady state (example)
HMO = y_all(:,1);
LfFe = y_all(:,2);

RBC_trace = y_all(:,13);
Hb_trace = Hb_estimator(RBC_trace);

figure;
plot(t_all, LfFe, 'LineWidth', 1.5);
xlabel('time (days)'); ylabel('HMO');
title('Model HMO (MATLAB translation)');
xlim([0 3])
ylim([0 0.004])
grid on;

figure;
plot(t_all, Hb_trace, 'LineWidth', 1.5);
xlabel('time (days)'); ylabel('Hb (g/dL)');
title('Model Hb trace (MATLAB translation)');
grid on;

figure;
plot(t_all, Hb_trace, 'LineWidth', 1.5);
xlabel('time (days)'); ylabel('Hb (g/dL)');
title('Model Hb trace (MATLAB translation)');
grid on;

figure;
semilogy(t_all, RBC_trace, 'LineWidth', 1.5);
xlabel('time (days)'); ylabel('RBC (log cells/kgbw)');
title('Model RBC trace (MATLAB translation)');
xlim([0 3])
grid on;

save('sim_standard_breastfeeding.mat','t_all','y_all','params','dose_HMO','dose_LfFe');
