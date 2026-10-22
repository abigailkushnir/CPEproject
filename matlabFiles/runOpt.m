clear;

parametersCPE;

ps0 = [P_NBCe2 P_NKKC1 P_NBCn1 P_NCBE P_AE2 P_K_a P_Cl_a]; % parameters being fitted

%% target concentrations and potentials
% targets: c_Na_i c_Cl_i c_K_i c_HCO_i Phi_i Phi_a
targets = [49.0 64.0 97.0 11.0 -0.0465 -0.0015];
targV_a = (targets(5) - targets(6))*F/(R*T);
targV_tj = -targets(6)*F/(R*T);
targX = targets(1) + targets(3) - targets(2) - targets(4); % from electroneutrality
weights = [1 1 1 1 1 1];

%% optimisation
opts = optimoptions('fmincon', 'Display', 'iter');
psFit = fmincon(@(ps) objFunc(ps, targets, weights), ps0, [], [], [], [], 0.0001*ps0, 10000*ps0, [], opts);

%% test the fitted parameters
syms c_Na_i c_Cl_i c_K_i c_HCO_i c_X_i V_a V_tj
vars = [c_Na_i c_Cl_i c_K_i c_HCO_i c_X_i V_a V_tj];
x0 = [targets(1) targets(2) targets(3) targets(4) targX targV_a targV_tj];
sol = vpasolve(modelEqs(vars, psFit) == 0, vars, x0);

if isempty(sol.c_Na_i)
    error('No solution for fitted ps');
end
x = double([sol.c_Na_i sol.c_Cl_i sol.c_K_i sol.c_HCO_i sol.c_X_i sol.V_a sol.V_tj]); % sym -> double

V_b = x(6) - x(7);
Phi_a = (V_b - x(6))*(R*T/F);
Phi_i = V_b*(R*T/F);

fprintf('Targets: c_Na_i = %.2f, c_Cl_i = %.2f, c_K_i = %.2f, c_HCO_i = %.2f, Phi_i = %.4f, Phi_a = %.4f\n', targets);
fprintf('Fitted:  c_Na_i = %.2f, c_Cl_i = %.2f, c_K_i = %.2f, c_HCO_i = %.2f, Phi_i = %.4f, Phi_a = %.4f\n', x(1:4), Phi_i, Phi_a);
fprintf('Fitted parameters: P_NBCe2 = %.4e, P_NKKC1 = %.4e, P_NBCn1 = %.4e, P_NCBE = %.4e, P_AE2 = %.4e, P_K_a = %.4e, P_Cl_a = %.4e\n', psFit);
