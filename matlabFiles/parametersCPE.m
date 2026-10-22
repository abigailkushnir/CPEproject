% apical concentrations (mM) (from ranges in Damkier et al. 2013)
c_Na_a = 152.0; 
c_Cl_a = 130.0; 
c_K_a = 3.0; 
c_HCO_a = 25.0;

% basal concentrations (mM) (from ranges in Damkier et al. 2013)
c_Na_b = 130.0; 
c_Cl_b = 114.0; 
c_K_b = 5.0; 
c_HCO_b = 21.0;

% apical and basal surface areas (m^2)
A_a = 1e-7; 
A_b = 1e-7;
A_tj = 1e-10;

% transporter intensities (mol/m^2/s)
P_ATP = 6e-6; 
P_NBCe2 = 6e-6; % being fitted 
P_NKKC1 = 1e-6; % being fitted 
P_NBCn1 = 1e-6; % being fitted 
P_NCBE = 1e-6;  % being fitted 
P_AE2 = 0.4e-6; % being fitted 

% channel permeabilities (m/s)
P_K_a = 9e-8; % being fitted
P_Cl_a = 7.4e-8; % being fitted 
P_Na_tj = 8e-5;
P_K_tj = 8e-5;
P_Cl_tj = 5e-5;
P_HCO_tj = 5e-5;

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
