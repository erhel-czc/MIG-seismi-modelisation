# Basic values
```python
#-------------------------------------------#
# Computational parameter definition
#-------------------------------------------#
dx=0.2                           # grid size
n=128                             # number of points
tol=1.0E-7                        # error tolerance
######################
# préférer tol=1.0E-12 pour les simulations, peut être baissé à 1.0E-7 pour des tests rapides
######################

nitrkmax=30                       # maximum number of iteration in a rkf step
nitmax=10000                        # maximum number of iterations 
hmin=1.0E-12                      # minimum time step
hmax=1.0E10                       # maximum time step (CFL for diffusion equation)
safe=0.8                         # safety factor for RKF iterations

#-------------------------------------#
# ND Mechanical parameter definition
#-------------------------------------#
alpha=0.2
eta=1.0E-8

#-------------------------------------------#
# Initial conditions (ND variables)
#-------------------------------------------#
v0=1.0       # initial normalized slip rate
th0=1/v0      # initial normalized state variable
t=0.0       # initial time
h=0.001     # initial time step
```