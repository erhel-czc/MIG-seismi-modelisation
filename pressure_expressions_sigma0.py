import math
import numpy as np


# Pressure expression from Journal of Geophysical Research: Solid Earth (P. Segall and S. Lu)
def P_article(t,pd,pnd):
    P0barre = (pd.Pinf/pd.sigma_n0)
    return P0barre * math.erfc(pnd.r / (2*math.sqrt(pnd.c * t)))

def dP_article(t,pd,pnd):
    P0barre = (pd.Pinf/pd.sigma_n0)
    return (P0barre / (2*t*math.sqrt(np.pi))) * (pnd.r / math.sqrt(pnd.c*t)) * np.exp(-(pnd.r/(2*math.sqrt(pnd.c*t)))**2)

# Linear pressure expression

def P_linear(t, pd, pnd):
    if t>=1e4:
        return pd.Pinf
    else :
        return (pd.Pinf/1e4)*t

def dP_linear(t, pd, pnd):
    if t>=1e4:
        return 0.0
    else :
        return pd.Pinf/1e4

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