import math
import numpy as np
from cycle_detection import find_cycle_start

taux = 0
path_without_pressure = "Results/PR_QDYN_RNS_modele_oriente/01"

# Pressure expression from Journal of Geophysical Research: Solid Earth (P. Segall and S. Lu)
def P_article(t,pd,pnd, delay = 0):
    if t < delay:
        return 0.0
    else :
        P0barre = (pd.Pinf/pd.sigma_n0)
        return P0barre * math.erfc(pnd.r / (2*math.sqrt(pnd.c * (t-delay))))

def dP_article(t,pd,pnd, delay = 0):
    P0barre = (pd.Pinf/pd.sigma_n0)
    if t < delay:
        return 0.0
    else :
        return (P0barre / (2*t*math.sqrt(np.pi))) * (pnd.r / math.sqrt(pnd.c*(t-delay)) * np.exp(-(pnd.r/(2*math.sqrt(pnd.c*(t-delay))))**2))

# Linear pressure expression

def P_linear(t, pd, pnd):
    t_max = 1e-3
    Pbarre = (pd.Pinf/pd.sigma_n0)
    if t>=t_max:
        return Pbarre
    else :
        return (Pbarre/t_max)*t+0.001

def dP_linear(t, pd, pnd):
    if t>=1e-3:
        return 0.001
    else :
        return pd.Pinf/(1e-3*pd.sigma_n0)

# Constant pressure

def P_constant(t, pd, pnd):
    if t<=2:
        return 0.0
    else:
        return pd.Pinf

def dP_constant(t, pd, pnd):
    return 0.0

# No pressure (P=0)
def P_none(t, pd, pnd, delay = find_cycle_start(path_without_pressure, taux)):
    return 0.0
def dP_none(t, pd, pnd, delay = find_cycle_start(path_without_pressure, taux)):
    return 0.0