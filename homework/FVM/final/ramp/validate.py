import numpy as np
from scipy.optimize import fsolve
import math

# Define parameters
theta_deg = 15  # given theta in degrees
theta = np.radians(theta_deg)  # convert to radians
M1 = 2.5  # given Mach number
gamma = 1.4  # specific heat ratio

# Define the theta-beta-M equation
def theta_beta_M(beta):
    beta_rad = np.radians(beta)
    numerator = M1**2 * (np.sin(beta_rad))**2 - 1
    denominator = M1**2 * (gamma + np.cos(2 * beta_rad)) + 2
    lhs = np.tan(theta)
    rhs = 2 * (1 / np.tan(beta_rad)) * (numerator / denominator)
    return lhs - rhs

# Initial guess for beta (in degrees)
beta_guess = 20  # reasonable guess for oblique shock

# Solve for beta
beta_solution = fsolve(theta_beta_M, beta_guess)


beta = np.radians(beta_solution)
print(f"Solution for beta: {beta_solution[0]:.4f} degrees")

# Pressure ratio p2/p1
p2_p1 = 1 + (2 * gamma / (gamma + 1)) * (M1**2 * math.sin(beta)**2 - 1)

# Density ratio rho2/rho1
rho2_rho1 = ((gamma + 1) * M1**2 * math.sin(beta)**2) / ((gamma - 1) * M1**2 * math.sin(beta)**2 + 2)

# Temperature ratio T2/T1
T2_T1 = (p2_p1) / (rho2_rho1)

# Post-shock Mach number M2
numerator = 1 + (gamma - 1) / 2 * M1**2 * math.sin(beta)**2
denominator = gamma * M1**2 * math.sin(beta)**2 - (gamma - 1) / 2
M2 = (1 / math.sin(beta - theta)) * math.sqrt(numerator / denominator)

p_in = 10e3
M_in = 2.5
T_in = 300
rho_in = p_in/287/T_in

print("density ps = ", rho2_rho1*rho_in)
print("temperature jump = ", T2_T1*T_in)
print("M2 = ", M2)