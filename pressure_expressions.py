import math
import numpy as np


# Pressure expression from Journal of Geophysical Research: Solid Earth (P. Segall and S. Lu)
def P_article(t,pd,pnd):
    P0barre = (pd.P0/pd.sigma_n)
    rbarre = (pd.r*pd.b_fric*pd.sigma_n) / (pd.mu*pd.dc)
    cbarre = pd.c * (pd.b_fric**2 * pd.sigma_n**2)/(pd.V_p * pd.mu**2 * pd.dc)
    return pd.sigma_n*P0barre * math.erfc(rbarre / (2*math.sqrt(cbarre * t)))

def dP_article(t,pd,pnd):
    P0barre = (pd.P0/pd.sigma_n)
    rbarre = (pd.r*pd.b_fric*pd.sigma_n) / (pd.mu*pd.dc)
    cbarre = pd.c * (pd.b_fric**2 * pd.sigma_n**2)/(pd.V_p * pd.mu**2 * pd.dc)
    return ((pd.V_p*pd.sigma_n*P0barre) / (pd.dc*2*t*math.sqrt(np.pi))) * (rbarre / math.sqrt(cbarre*t)) * np.exp(-(rbarre/(2*math.sqrt(cbarre*t)))**2)

# Linear pressure expression

def P_linear(t, pd, pnd):
    if t>=1e4:
        return pd.P0
    else :
        return (pd.P0/1e4)*t

def dP_linear(t, pd, pnd):
    if t>=1e4:
        return 0.0
    else :
        return pd.P0/1e4

# Constant pressure

def P_constant(t, pd, pnd):
    return pd.P0

def dP_constant(t, pd, pnd):
    return 0.0

# No pressure (P=0)
def P_none(t, pd, pnd):
    return 0.0
def dP_none(t, pd, pnd):
    return 0.0