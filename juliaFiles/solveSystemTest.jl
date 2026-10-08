using Pkg

Pkg.activate("./")
# Pkg.add("NonlinearSolve")
# Pkg.add("Plots")
# Pkg.add("ModelingToolkit")
using NonlinearSolve
using ModelingToolkit

@variables Phi_i Phi_a Phi_b V_a V_b # interior, apical, and basal potentials and apical and basal potential differences
@variables c_Na_i c_Cl_i c_K_i c_HCO_i c_X_i # interior concentrations of Sodium, Cloride, Potassium, Bicorbonate, and large proteins, respectively

@parameters c_Na_a c_Cl_a c_K_a c_HCO_a c_Na_b c_Cl_b c_K_b c_HCO_b # apical and basal concentrations
@parameters A_a A_b # apical and basal surface areas (m^2)
@variables K_Na K_K # dissaociation constants for Sodium and Potassium (mM)
@parameters P_ATP P_NBCe2 P_NKKC1 P_NBCn1 P_NCBE P_AE2 # transporter intensities (mol/m^2/s) 
@parameters P_K P_Cl # channel permeabilities (m/s)
@parameters σ_Na σ_K σ_Cl σ_HCO σ_X # reflection coefficients for Sodium, Potassium, Chloride, Bicarbonate, and large proteins, respectively
@parameters p_i p_b p_a # interior, basal, and apical pressures
@parameters R T F # gas constant, temperature, and Faraday's constant
@parameters K_a K_b # apical and basal hydraulic conductivities (m/s/Pa)

# Apical Fluxes #
J_ATP = A_a * P_ATP * (c_Na_i/(K_Na + c_Na_i))^3 * (c_K_a/(K_K + c_K_a))^2 # ATPase flux
J_NBCe2 = A_a * P_NBCe2 * (log((c_Na_i * c_HCO_i^3)/(c_Na_a * c_HCO_a^3)) - 2*V_a) # NBCe2 flux
J_NKCC1 = A_a * P_NKKC1 * log((c_Na_i * c_K_i * c_Cl_i^2)/(c_Na_a * c_K_a * c_Cl_a^2)) # NKCC1 flux
J_K = A_a * P_K * V_a * ((c_K_i - c_K_a * exp(-V_a))/(1 - exp(-V_a))) # Potassium channel flux
J_Cl = A_a * P_Cl * (-1) * V_a * ((c_Cl_i - c_Cl_a * exp(V_a))/(1 - exp(V_a))) # Chloride channel flux
# Basal Fluxes #
J_NBCn1 = A_b * P_NBCn1 * log((c_Na_b * c_HCO_b)/(c_Na_i * c_HCO_i)) # NBCn1 flux
J_NCBE = A_b * P_NCBE * log((c_Na_b * c_Cl_i * c_HCO_b^2)/(c_Na_i * c_Cl_b * c_HCO_i^2)) # NCBE flux
J_AE2 = A_b * P_AE2 * log((c_Cl_b * c_HCO_i)/(c_Cl_i * c_HCO_b)) # AE2 flux
# Water Fluxes #
Q_a = A_a * K_a *((p_i - p_a) - R*T*(σ_Na*(c_Na_i - c_Na_a) + σ_K*(c_K_i - c_K_a) + σ_Cl*(c_Cl_i - c_Cl_a) + σ_HCO*(c_HCO_i - c_HCO_a) + σ_X*(c_X_i))) # apical water flux
Q_b = A_b * K_b *((p_i - p_b) - R*T*(σ_Na*(c_Na_i - c_Na_b) + σ_K*(c_K_i - c_K_b) + σ_Cl*(c_Cl_i - c_Cl_b) + σ_HCO*(c_HCO_i - c_HCO_b) + σ_X*(c_X_i))) # basal water flux

# Ion Fluxes #
J_Na_a = 3 * J_ATP + J_NBCe2 + J_NKCC1 # apical Sodium flux
J_Na_b = - J_NBCn1 - J_NCBE # basal Sodium flux
J_Cl_a = 2 * J_NKCC1 + J_Cl # apical Chloride flux
J_Cl_b = J_NCBE - J_AE2 # basal Chloride flux
J_K_a = -2 * J_ATP + J_NKCC1 + J_K # apical Potassium flux
J_K_b = 0 # basal Potassium flux
J_HCO_a = 3 * J_NBCe2 # apical Bicarbonate flux
J_HCO_b = - J_NBCn1 - 2 * J_NCBE + J_AE2 # basal Bicarbonate flux

# Balance Equations #
eqs = [
    V_a ~ (F/(R*T))*(Phi_i - Phi_a) # apical potential difference
    V_b ~ (F/(R*T))*(Phi_i - Phi_b) # basal potential difference
    Phi_a ~ 0 # apical potential
    Phi_b ~ 0 # basal potential
    K_Na ~ 0.2 * (1 + c_K_i/8.33) # Sodium dissociation constant
    K_K ~ 0.1 * (1 + c_Na_a/18.5) # Potassium dissociation constant

    0 ~ (J_Na_a + J_Na_b) / (A_a * 1e-6) # Sodium balance
    0 ~ (J_Cl_a + J_Cl_b) / (A_a * 1e-6) # Chloride balance
    0 ~ (J_K_a + J_K_b) / (A_a * 1e-6) # Potassium balance
    0 ~ (J_HCO_a + J_HCO_b) / (A_a * 1e-6) # Bicarbonate balance
    0 ~ (Q_a + Q_b) / (A_a * K_a * R * T) # Water balance
    0 ~ c_Na_i + c_K_i - c_Cl_i - c_HCO_i - c_X_i # electroneutrality

]

@mtkcompile ns = System(eqs)

ps = [
# Set parameter values #
    # Apical concentrations #
c_Na_a => 130.0,
c_Cl_a => 114.0,
c_K_a => 5.0,
c_HCO_a => 21.0,
    # Basal concentrations #
c_Na_b => 130.0,
c_Cl_b => 114.0,
c_K_b => 5.0,
c_HCO_b => 21.0,
    # Apical and basal surface areas #
A_a => 6 * 10^(-4),
A_b => 6 * 10^(-4),
    # Transporter inensities #
P_ATP => 6 * 10^(-6),
P_NBCe2 => 6 * 10^(-6),
P_NKKC1 => 1 * 10^(-6),
P_NBCn1 => 1 * 10^(-6),
P_NCBE => 1 * 10^(-6),
P_AE2 => 0.4 * 10^(-6),
    # Channel permeabilities #
P_K => 5 * 10^(-8),
P_Cl => 6 * 10^(-8),
    # Reflection coefficients #
σ_Na => 1,
σ_K => 1,
σ_Cl => 1,
σ_HCO => 1,
σ_X => 1,
    # Pressures #
p_i => 0,
p_a => 0,
p_b => 0,
    # Gas constant, temperature, and Faraday's constant #
R => 8.314,
T => 310,
F => 96485.3321,
    # Apical and basal hydraulic conductivities #
K_a => 2 * 10^(-10),
K_b => 2 * 10^(-10),
]
guesses = [c_Na_i => 1.0, c_Cl_i => 1.0, c_K_i => 1.0, c_HCO_i => 1.0, c_X_i => 1.0, Phi_i => 1, Phi_a => 0.0, Phi_b => 0.0, V_a => -0.1, V_b => -0.1, K_Na => 0.1, K_K => 0.1]

prob = NonlinearProblem(ns, vcat(guesses, ps))
sol  = solve(prob, TrustRegion())
println(sol.retcode)

