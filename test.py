import numpy as np
import matplotlib as plt

# velocity in m/s
V = 1 
# diameter in m
D = 0.01 
# temperature in K
T = 300 

# liquid water properties at 1 atm
def rho(T):
    return 1000 * (1 - ((T + 288.9414)/(508929.2*(T + 68.12963))) * (T - 3.9863)**2)
def mu(T):
    return 2.414e-5 * 10**(247.8/(T + 133.15))
def Cp(T):
    4179 + 0.1*(T-20)
def k(T):
    4179 + 0.1*(T-20)
def Re(T):   
    return rho(T)*V*D/mu(T)
def Pr(T):
    return mu(T)*Cp(T)/k(T)

# Correlations 
# je fais genre j'ai modifié
# Re > 10000 and 0.7 < Pr < 160 and L/D > 10
def dittus_bolter(Re,Pr):
    if Re < 8000:
        raise ValueError("Reynolds number must be greater than 10000")
    else:
        alpha = 0.8
        beta = 0.3
        Nu = 0.023 * Re**alpha * Pr**beta
    return Nu

# 3000 < Re < 5 000 000 and 0.5 < Pr < 2000
def petukhov(Re):
    if Re < 2300:
        f = 64/Re
    elif 2300 < Re < 4000:
        f = (1 - ((Re - 2300) / (4000 - 2300))) * 64 / Re + ((Re - 2300) / (4000 - 2300)) * (0.79 * np.log(Re) - 1.64)**(-2)
    else:
        f = (0.79 * np.log(Re) - 1.64)**(-2)
    return f

def colebrook(Re, eps, D):
    if Re < 2300:
        f = 64/Re
    elif Re < 4000:
        f_lam = 64/Re
        f_turb = 0.25/(np.log10(eps/(3.7*D) + 5.74/Re**0.9)**2)
        w = (Re-2300)/(4000-2300)
        f = (1-w)*f_lam + w*f_turb  # transition (approx)
    else:
        f = 0.25/(np.log10(eps/(3.7*D) + 5.74/Re**0.9)**2)  # Swamee-Jain
    return f

def gnielinski(Re,Pr):
    if Re < 3000:
        Nu = 4.36
    else:
        f = petukhov(Re)
        Nu = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * (f / 8)**0.5 * (Pr**(2/3) - 1))
    return Nu

# If the spread between the surface and the bulk temperature is high (50..)
# The Sieder-Tate correlation can be used to account for the variation of viscosity with temperature.


# Re > 10000 and 0.7 < Pr < 16700 and L/D > 10
def Sieder_Tate(Re,Pr,mu_surf,mu_bulk):
    if Re < 10000:
        raise ValueError("Reynolds number must be greater than 10000")
    else:
        alpha = 4/5
        beta = 1/3
        gamma = 0.14
        Nu = 0.027*Re**alpha*Pr**beta*(mu_bulk/mu_surf)**gamma
    return Nu

# Exchange coefficients for different geometries can be calculated using the appropriate correlations.
def h(T):
    Nu_dittus, Nu_petukhov, Nu_gnielinski
    k_val = k(T)
    h_dittus = Nu_dittus * k_val / D
    h_petukhov = Nu_petukhov * k_val / D
    h_gnielinski = Nu_gnielinski * k_val / D
    return h_dittus, h_petukhov, h_gnielinski

# Class (straight pipes, annuli pipes, elbow pipes, etc)

#Straight, define Q if you want to give power to the fluide. Or if you want to extract. 
L_pipe = 1 # length of the pipe in m
rugosity = 0.0001 # pipe roughness in m

def straight_pipe(T,Q,L,f):
    h_val = h(T)
    Tin = T
    Ts = Q / (h_val * np.pi * D) + Tin
    Tout = Tin + Q / (rho(T)*Cp(T)*V*np.pi*D**2/4)
    dp = f * L / D * rho(T) * V**2 / 2
    return Ts, Tout, dp

def elbow_pipe(T,K):
    rho_val = rho(T)
    dp_elbow = K * rho_val * V**2 / 2
    return dp_elbow

def valve_pipe(T,K2):
    rho_val = rho(T)
    dp_valve = K2 * rho_val * V**2 / 2
    return dp_valve
