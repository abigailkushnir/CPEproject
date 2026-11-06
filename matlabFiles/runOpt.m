clear;

parametersCPE;

ps0 = [P_NBCe2 P_NKKC1 P_NBCn1 P_NCBE P_AE2 P_K_a P_Cl_a P_K_b P_Na_b P_K_tj P_Cl_tj P_Na_tj P_HCO_tj]; % parameters being fitted

%% target concentrations and potentials
% targets: c_Na_i c_Cl_i c_K_i c_HCO_i Phi_i Phi_a
targets = [49.0 64.0 97.0 11.0 -0.0465 -0.0015];
targV_a = (targets(5) - targets(6))*F/(R*T);
targV_tj = -targets(6)*F/(R*T);
targX = targets(1) + targets(3) - targets(2) - targets(4); % from electroneutrality
weights = [1 1 1 1 1 1];

%% optimisation
lb = 0.1*ps0; % lower bounds
ub = 10*ps0; % upper bounds
opts = optimoptions('fmincon', 'Display', 'iter'); %, 'StepTolerance', 1e-16, 'OptimalityTolerance', 1e-8);
% psFit = fmincon(@(ps) objFunc(ps, targets, weights), ps0, [], [], [], [], lb, ub, [], opts);
qFit = fmincon(@(q) objFunc(10.^q, targets, weights), log10(ps0), [], [], [], [], log10(lb), log10(ub), [], opts);
psFit = 10.^qFit;

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

% use fitted parameters in expressions
P_NBCe2 = psFit(1); P_NKKC1 = psFit(2); P_NBCn1 = psFit(3); P_NCBE = psFit(4); P_AE2 = psFit(5); P_K_a = psFit(6); P_Cl_a = psFit(7); P_K_b = psFit(8); P_Na_b = psFit(9); P_K_tj = psFit(10); P_Cl_tj = psFit(11); P_Na_tj = psFit(12); P_HCO_tj = psFit(13);
expressions;

fprintf('Targets: c_Na_i = %.2f, c_Cl_i = %.2f, c_K_i = %.2f, c_HCO_i = %.2f, Phi_i = %.4f, Phi_a = %.4f\n', targets);
fprintf('Fitted:  c_Na_i = %.2f, c_Cl_i = %.2f, c_K_i = %.2f, c_HCO_i = %.2f, Phi_i = %.4f, Phi_a = %.4f\n', x(1:4), Phi_i, Phi_a);
fprintf('Fitted parameters: \n \tChannels: P_NBCe2 = %.4e, P_NKKC1 = %.4e, P_NBCn1 = %.4e, P_NCBE = %.4e, P_AE2 = %.4e, \n \tPermeabilities: P_K_a = %.4e, P_Cl_a = %.4e\n, P_K_b = %.4e, P_Na_b = %.4e, P_K_tj = %.4e, P_Cl_tj = %.4e, P_Na_tj = %.4e, P_HCO_tj = %.4e', psFit);
fprintf('\nI_a: %.4e, I_tj: %.4e\n', double(subs(I_a, vars, x)), double(subs(I_tj, vars, x)));
fprintf('\nJ_K_channel_b: %.4e, J_Na_channel_b: %.4e\n', double(subs(J_K_channel_b, vars, x)), double(subs(J_Na_channel_b, vars, x)));