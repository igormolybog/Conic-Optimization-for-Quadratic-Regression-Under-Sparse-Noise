% mpc = case1354pegase;
mpc = case9;
len = 60;
rng(42)

active_process = path_generation(mpc.bus(:,3), len);
reactive_process = path_generation(mpc.bus(:,4), len);

V = zeros(size(mpc.bus, 1), len+1);
V_mp = zeros(size(mpc.bus, 1), len+1);
distance = zeros(1, len+1);
[vV(:, 1), vV_mp(:, 1), vdistance(:, 1)] = socp_experiment(mpc, 5e-1);
% [yV(:, 1), yV_mp(:, 1), ydistance(:, 1)] = socp_experiment(mpc, 5e-1);
vV(:, 1) = vV_mp(:, 1);

for i = 1:len
    mpc.bus(:,3) = active_process(:, i);
    mpc.bus(:,4) = reactive_process(:, i);
    [vV(:, i+1), vV_mp(:, i+1), vdistance(:, i+1)] = socp_experiment(mpc, 5e-3, vV(:, i));
%     [yV(:, i+1), yV_mp(:, i+1), ydistance(:, i+1)] = socp_experiment(mpc, 5e-1);
end





% 
% res = runpf(mpc);
% VM = 8; % MATPOWER 7.0 p. 149 https://matpower.org/docs/manual.pdf
% VA = 9;
% voltageMagnitudes = res.bus(:, VM);
% voltageAngles = res.bus(:, VA);