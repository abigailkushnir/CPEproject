% dissociation constants
K_Na = 0.2*(1 + c_K_i/8.33); K_K = 0.1*(1 + c_Na_a/18.5);

% apical
J_ATP = A_a*P_ATP*(c_Na_i/(K_Na + c_Na_i))^3*(c_K_a/(K_K + c_K_a))^2; % ATPase
J_NBCe2 = A_a*P_NBCe2*(log((c_Na_i*c_HCO_i^3)/(c_Na_a*c_HCO_a^3)) - 2*V_a); % NBCe2
J_NKCC1 = A_a*P_NKKC1*log((c_Na_i*c_K_i*c_Cl_i^2)/(c_Na_a*c_K_a*c_Cl_a^2)); % NKCC1
J_K_channel_a = A_a*P_K_a*V_a*(c_K_i - c_K_a*exp(-V_a))/(1 - exp(-V_a)); % K channel
J_Cl_channel_a = -A_a*P_Cl_a*V_a*(c_Cl_i - c_Cl_a*exp(V_a))/(1 - exp(V_a)); % Cl channel

% basal
J_NBCn1 = A_b*P_NBCn1*log((c_Na_b*c_HCO_b)/(c_Na_i*c_HCO_i)); % NBCn1
J_NCBE = A_b*P_NCBE*log((c_Na_b*c_Cl_i*c_HCO_b^2)/(c_Na_i*c_Cl_b*c_HCO_i^2)); % NCBE
J_AE2 = A_b*P_AE2*log((c_Cl_b*c_HCO_i)/(c_Cl_i*c_HCO_b)); % AE2
J_K_channel_b = A_b*P_K_b*(V_a - V_tj)*(c_K_i - c_K_b*exp(-(V_a - V_tj)))/(1 - exp(-(V_a - V_tj))); % K channel
J_Na_channel_b = A_b*P_Na_b*(V_a - V_tj)*(c_Na_i - c_Na_b*exp(-(V_a - V_tj)))/(1 - exp(-(V_a - V_tj))); % Na channel

% water
Q_a = A_a*K_a*((p_i - p_a) - R*T*(sigma_Na*(c_Na_i - c_Na_a) + sigma_K*(c_K_i - c_K_a) + sigma_Cl*(c_Cl_i - c_Cl_a) + sigma_HCO*(c_HCO_i - c_HCO_a) + sigma_X*c_X_i));
Q_b = A_b*K_b*((p_i - p_b) - R*T*(sigma_Na*(c_Na_i - c_Na_b) + sigma_K*(c_K_i - c_K_b) + sigma_Cl*(c_Cl_i - c_Cl_b) + sigma_HCO*(c_HCO_i - c_HCO_b) + sigma_X*c_X_i));

% ion fluxes
J_Na_a = 3*J_ATP + J_NBCe2 + J_NKCC1;
J_Na_b = -J_NBCn1 - J_NCBE + J_Na_channel_b;
J_Na_tj = A_tj*P_Na_tj*V_tj*((c_Na_b - c_Na_a*exp(-V_tj))/(1 - exp(-V_tj)));

J_Cl_a = 2*J_NKCC1 + J_Cl_channel_a;
J_Cl_b = J_NCBE - J_AE2;
J_Cl_tj = -A_tj*P_Cl_tj*V_tj*((c_Cl_b - c_Cl_a*exp(V_tj))/(1 - exp(V_tj)));

J_K_a = -2*J_ATP + J_NKCC1 + J_K_channel_a;
J_K_b = J_K_channel_b;
J_K_tj = A_tj*P_K_tj*V_tj*((c_K_b - c_K_a*exp(-V_tj))/(1 - exp(-V_tj)));

J_HCO_a = 3*J_NBCe2;
J_HCO_b = -J_NBCn1 - 2*J_NCBE + J_AE2;
J_HCO_tj = -A_tj*P_HCO_tj*V_tj*((c_HCO_b - c_HCO_a*exp(V_tj))/(1 - exp(V_tj)));

% current equations
I_a = F*(J_ATP - 2*J_NBCe2 + J_K_channel_a - J_Cl_channel_a);
% I_b = 0;
I_tj = F*(J_Na_tj + J_K_tj - J_Cl_tj - J_HCO_tj);