# Valeurs de référence des paramètres dimensionnels utilisés dans les simulations

```python
shear=2.0E10    # shear modulus (Pa) 
rho_roc=2700.0  # rock density (kg/m3)
lenght_fault=1.0E3  # fault lenght (m)
depth_fault=3.0E3   # fault depth (m)
a_fric=0.005      # direct effect coefficient
b_fric=0.01       # evolution effect coefficient
dc=1.0E-4           # critical slip distance (m)
V_p=1.0E-9        # tectonic speed (m/s)


#-------------------------------------------#
# Initial conditions (ND variables)
#-------------------------------------------#
v=2       # initial normalized slip rate
th=1/v      # initial normalized state variable
sigma_n=1.0  # initial normal stress (ND)

t=0.0       # initial time
h=0.001     # initial time step

psi=0.7
f0=0.6
```

## Best values
```
shear=2.0E10    # shear modulus (Pa)
rho_roc=2700.0  # rock density (kg/m3)
lenght_fault=1.0E3  # fault lenght (m)
depth_fault=3.0E3   # fault depth (m)
a_fric=0.005    # direct effect coefficient
b_fric=0.01       # evolution effect coefficient
dc=1.0E-3           # critical slip distance (m)
V_p=1.0E-7        # tectonic speed (m/s)
r_real = 1.0e2 # distance to the injection point (m)
c_real = 6.8e-4  # hydraulic diffusivity (m2/s)
Pinf = 2.5e8     # injection pressure (Pa)
```