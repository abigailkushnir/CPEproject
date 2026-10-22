parametersCPE;

syms c_Na_i c_Cl_i c_K_i c_HCO_i c_X_i V_a V_tj

K_Na = 0.2*(1 + c_K_i/8.33); K_K = 0.1*(1 + c_Na_a/18.5); % dissociation constants

% apical
J_ATP = A_a*P_ATP*(c_Na_i/(K_Na + c_Na_i))^3*(c_K_a/(K_K + c_K_a))^2; % ATPase
J_NBCe2 = A_a*P_NBCe2*(log((c_Na_i*c_HCO_i^3)/(c_Na_a*c_HCO_a^3)) - 2*V_a); % NBCe2
J_NKCC1 = A_a*P_NKKC1*log((c_Na_i*c_K_i*c_Cl_i^2)/(c_Na_a*c_K_a*c_Cl_a^2)); % NKCC1
J_K = A_a*P_K_a*V_a*(c_K_i - c_K_a*exp(-V_a))/(1 - exp(-V_a)); % K channel
J_Cl = -A_a*P_Cl_a*V_a*(c_Cl_i - c_Cl_a*exp(V_a))/(1 - exp(V_a)); % Cl channel

% basal
J_NBCn1 = A_b*P_NBCn1*log((c_Na_b*c_HCO_b)/(c_Na_i*c_HCO_i)); % NBCn1
J_NCBE = A_b*P_NCBE*log((c_Na_b*c_Cl_i*c_HCO_b^2)/(c_Na_i*c_Cl_b*c_HCO_i^2)); % NCBE
J_AE2 = A_b*P_AE2*log((c_Cl_b*c_HCO_i)/(c_Cl_i*c_HCO_b)); % AE2

% water
Q_a = A_a*K_a*((p_i - p_a) - R*T*(sigma_Na*(c_Na_i - c_Na_a) + sigma_K*(c_K_i - c_K_a) + sigma_Cl*(c_Cl_i - c_Cl_a) + sigma_HCO*(c_HCO_i - c_HCO_a) + sigma_X*c_X_i));
Q_b = A_b*K_b*((p_i - p_b) - R*T*(sigma_Na*(c_Na_i - c_Na_b) + sigma_K*(c_K_i - c_K_b) + sigma_Cl*(c_Cl_i - c_Cl_b) + sigma_HCO*(c_HCO_i - c_HCO_b) + sigma_X*c_X_i));

% ion fluxes
J_Na_a = 3*J_ATP + J_NBCe2 + J_NKCC1;  
J_Na_b = -J_NBCn1 - J_NCBE;
J_Na_tj = A_tj*P_Na_tj*V_tj*((c_Na_b - c_Na_a * exp(-V_tj))/(1 - exp(-V_tj)));

J_Cl_a = 2*J_NKCC1 + J_Cl;             
J_Cl_b = J_NCBE - J_AE2;
J_Cl_tj = -A_tj*P_Cl_tj*V_tj*((c_Cl_b - c_Cl_a*exp(V_tj))/(1 - exp(V_tj)));

J_K_a = -2*J_ATP + J_NKCC1 + J_K;      
J_K_b = 0;
J_K_tj = A_tj*P_K_tj*V_tj*((c_K_b - c_K_a*exp(-V_tj))/(1 - exp(-V_tj)));

J_HCO_a = 3*J_NBCe2;                   
J_HCO_b = -J_NBCn1 - 2*J_NCBE + J_AE2;
J_HCO_tj = -A_tj*P_HCO_tj*V_tj*((c_HCO_b - c_HCO_a*exp(V_tj))/(1 - exp(V_tj)));

% current equations
I_a = F*(J_ATP - 2*J_NBCe2 + J_K - J_Cl);
I_tj = F*(J_Na_tj + J_K_tj - J_Cl_tj - J_HCO_tj);

% scales
J_scale = A_a*P_ATP; % ion flux scale (mol/s)
Q_scale = A_a*K_a*R*T; % water flux per mM osmotic difference (m^3/s)

% balance equations
eq1 = (J_Na_a + J_Na_b)/J_scale == 0; % sodium
eq2 = (J_Cl_a + J_Cl_b)/J_scale == 0; % chloride
eq3 = (J_K_a + J_K_b)/J_scale == 0; % potassium
eq4 = (J_HCO_a + J_HCO_b)/J_scale == 0; % bicarbonate
eq5 = (Q_a + Q_b)/Q_scale == 0; % water
eq6 = c_Na_i + c_K_i - c_Cl_i - c_HCO_i - c_X_i == 0; % electroneutrality
eq7 = (I_a + I_tj)/(F*J_scale) == 0; % open circuit

eqns = [eq1 eq2 eq3 eq4 eq5 eq6 eq7];
vars = [c_Na_i c_Cl_i c_K_i c_HCO_i c_X_i V_a V_tj];
incon = [10, 40, 130, 20, 80, -1.5, -0.05]; % V_tj starts > 0, away from the 0/0 at V_tj = 0
clear sol;

sol = vpasolve(eqns, vars, incon);

if isempty(sol.c_Na_i)==0 && isreal(sol.c_Na_i)==1 && isreal(sol.V_a)==1

c_Na_i = double(sol.c_Na_i);
c_Cl_i = double(sol.c_Cl_i);
c_K_i = double(sol.c_K_i);
c_HCO_i = double(sol.c_HCO_i);
c_X_i = double(sol.c_X_i);
V_a = double(sol.V_a);
V_tj = double(sol.V_tj);
K_Na = 0.2*(1 + c_K_i/8.33);
K_K = 0.1*(1 + c_Na_a/18.5);

fprintf('Na: %.2f\n',c_Na_i);
fprintf('Cl: %.2f\n',c_Cl_i);
fprintf('K: %.2f\n', c_K_i);
fprintf('HCO3: %.2f\n', c_HCO_i);
fprintf('X: %.2f\n', c_X_i);
fprintf('V_a: %.4f\n', V_a);
fprintf('V_tj: %.4f\n', V_tj);
V_b = V_a - V_tj;
fprintf('Phi_a*1000: %.4f\n', (V_b - V_a)*(R*T/F)*1000);
fprintf('Phi_i*1000: %.4f\n', (V_b)*(R*T/F)*1000);

else
   fprintf("failed\n")
   c_Na_i = NaN; c_Cl_i = NaN; c_K_i = NaN; c_HCO_i = NaN; c_X_i = NaN;
   V_a = NaN; V_b = NaN; Phi_i = NaN; Phi_a = NaN; Phi_b = NaN; K_Na = NaN; K_K = NaN;
end
