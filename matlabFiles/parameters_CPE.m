% apical concentrations (mM)
c_Na_a = 130.0; 
c_Cl_a = 114.0; 
c_K_a = 5.0; 
c_HCO_a = 21.0;

% basal concentrations (mM)
c_Na_b = 130.0; 
c_Cl_b = 114.0; 
c_K_b = 5.0; 
c_HCO_b = 21.0;

% apical and basal surface areas (m^2)
A_a = 6e-4; 
A_b = 6e-4;
A_tj = 2e-4;

% transporter intensities (mol/m^2/s)
P_ATP = 6e-6; 
P_NBCe2 = 6e-6; 
P_NKKC1 = 1e-6;
P_NBCn1 = 1e-6; 
P_NCBE = 1e-6; 
P_AE2 = 0.4e-6;

% channel permeabilities (m/s)
P_K_a = 5e-8; 
P_Cl_a = 6e-8;
P_Na_tj = 5e-8;
P_K_tj = 5e-8;
P_Cl_tj = 5e-8;
P_HCO_tj = 5e-8;

% reflection coefficients
sigma_Na = 1; 
sigma_K = 1; 
sigma_Cl = 1; 
sigma_HCO = 1; 
sigma_X = 1;

% interior, apical, and basal pressures
p_i = 0; 
p_a = 0; 
p_b = 0;

% gas constant, temperature, and Faraday's constant
R = 8.314; 
T = 310; 
F = 96485.3321;

% apical and basal hydraulic conductivities (m/s/Pa)
K_a = 2e-10; 
K_b = 2e-10;
