function y = modelEqs(vars, ps)

    parametersCPE;

    % parameters being varied
    P_NBCe2 = ps(1); P_NKKC1 = ps(2); P_NBCn1 = ps(3); P_NCBE = ps(4); P_AE2 = ps(5); P_K_a = ps(6); P_Cl_a = ps(7); P_K_b = ps(8); P_Na_b = ps(9); P_K_tj = ps(10); P_Cl_tj = ps(11); P_Na_tj = ps(12); P_HCO_tj = ps(13);
    c_Na_i = vars(1); c_Cl_i = vars(2); c_K_i = vars(3); c_HCO_i = vars(4); c_X_i = vars(5); V_a = vars(6); V_tj = vars(7);

    expressions;

    % scales
    J_scale = A_a*P_ATP; % ion flux scale (mol/s)
    Q_scale = A_a*K_a*R*T; % water flux per mM osmotic difference (m^3/s)

    % balance equations
    y(1) = (J_Na_a + J_Na_b)/J_scale;                  % sodium
    y(2) = (J_Cl_a + J_Cl_b)/J_scale;                  % chloride
    y(3) = (J_K_a + J_K_b)/J_scale;                    % potassium
    y(4) = (J_HCO_a + J_HCO_b)/J_scale;                % bicarbonate
    y(5) = (Q_a + Q_b)/Q_scale;                        % water
    y(6) = c_Na_i + c_K_i - c_Cl_i - c_HCO_i - c_X_i;  % electroneutrality
    y(7) = (I_a + I_tj)/(F*J_scale);                  % open circuit

end