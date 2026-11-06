function L = objFuncFsolve(ps, targets, weights)

    parametersCPE;

    % targets = [c_Na_i c_Cl_i c_K_i c_HCO_i Phi_i Phi_a]
    targV_a = (targets(5) - targets(6))*F/(R*T);
    targV_tj = -targets(6)*F/(R*T);
    targX = targets(1) + targets(3) - targets(2) - targets(4); % c_X_i from electroneutrality

    x0 = [targets(1) targets(2) targets(3) targets(4) targX targV_a targV_tj]; % initial guess near targets

    % numeric solve (much faster than vpasolve)
    solveOpts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);
    try
        [x, ~, exitflag] = fsolve(@(v) modelEqs(v, ps), x0, solveOpts);
    catch
        exitflag = -1; % e.g. complex values from log of negative argument
    end
    if exitflag <= 0
        L = 1e10; % no solution found
        return
    end

    V_b = x(6) - x(7);
    Phi_a = (V_b - x(6))*(R*T/F);
    Phi_i = V_b*(R*T/F);

    r = [x(1:4) Phi_i Phi_a] - targets;

    if ~isreal(x) || any(x(1:5) <= 0)
        L = 1e10; % complex or neg concentrations
    else
        L = sum(weights.*((r)./targets).^2); % fmincon finding local min better when rel. error is used
        % L = sum(weights.*(r).^2);
    end

end
