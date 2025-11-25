import numpy as np
from result import Result

#####################################
# Parameters
#####################################

#-------------------------------------#
# Dimensional Mechanical parameter definition
#-------------------------------------#
k_rigidity=3.0E7 # rigidity (Pa)
a_fric=0.005      # direct effect coefficient
b_fric=0.01       # evolution effect coefficient
eta_visc=1.0E6   # viscosity (Pa.s)
sigma_n0=1.0E8    # normal stress (Pa)
dc=1.0E-4           # critical slip distance (m)
V_p=1.0E-9        # tectonic speed (m/s)

#-------------------------------------#
# ND Mechanical parameter definition
#-------------------------------------#
a=0.6
eta=1.0E-11
k=0.41

#-------------------------------------------#
# Computational parameter definition
#-------------------------------------------#
tol=1.0E-10   # error tolerance
nitrkmax=30
nitmax=10000  # maximum number of iterations
hmin=1.0E-12  # minimum time step
hmax=1.0E10   # maximum time step
safe=0.8      # safety factor for RKF iterations

#-------------------------------------------#
# Initial conditions (ND variables)
#-------------------------------------------#
v=0.1       # initial normalized slip rate
th=1/v      # initial normalized state variable
sigma_n=1.0  # initial normal stress (ND)

t=0.0       # initial time
h=0.001     # initial time step

psi=0.7
f0=0.6

phi=np.log(v)
nu=np.log(th)

class ParamMec:
    "Dimensional Mechanical parameters"

    def __init__(self, shear, rho_rock, lenght_fault, depth_fault, a_fric, b_fric, dc, V_p):
        self.shear=shear
        self.rho_rock=rho_rock
        self.lenght_fault=lenght_fault
        self.depth_fault=depth_fault
        self.a_fric=a_fric
        self.b_fric=b_fric        
        self.dc=dc
        self.V_p=V_p
        self.sigma_n0=rho_rock*9.81*depth_fault  # lithospheric stress
        self.k_rigidity= shear/lenght_fault  # rigidity
        self.eta_visc= np.sqrt(shear*rho_rock)/2. # viscosity

        

class NdParamMec:
    "Non dimensional Mechanical parameters"

    def __init__(self, a, eta, k, psi, f0, b):
        self.a=a
        self.eta=eta
        self.k=k
        self.psi=psi
        self.f0=f0
        self.b=b

class ParamComp:
    "Computational parameters"

    def __init__(self, tol, nitrkmax, nitmax, hmin, hmax, safe):
        self.tol=tol # error tolerance
        self.nitrkmax=nitrkmax # maximum number of iteration in a rkf step
        self.nitmax=nitmax # maximum number of iterations
        self.hmin=hmin # minimum time step
        self.hmax=hmax # maximum time step (CFL for diffusion equation)
        self.safe=safe




#-------------------------------------#
# Dimensional Mechanical parameter definition
#-------------------------------------#

pd = ParamMec(k_rigidity=k_rigidity, a_fric=a_fric, b_fric=b_fric, eta_visc=eta_visc, sigma_n0=sigma_n0, dc=dc, V_p=V_p)

#-------------------------------------#
# ND Mechanical parameter definition
#-------------------------------------#

pnd=NdParamMec(a = pd.a_fric/pd.b_fric, k = pd.k_rigidity*pd.dc/(pd.sigma_n0*pd.b_fric), eta = pd.eta_visc*pd.V_p/(pd.b_fric*pd.sigma_n0), psi=psi, f0=f0, b=pd.b_fric)

f=pnd.f0 + pnd.a*pnd.b*np.log(v) + pnd.b*np.log(th)  # initial frictional resistance (ND)

cpsi=np.cos(pnd.psi)
spsi=np.sin(pnd.psi)

#-------------------------------------------#
# Computational parameter definition
#-------------------------------------------#
pc = ParamComp(tol, nitrkmax, nitmax, hmin, hmax, safe)



def f_rns(phi, nu, pnd):

    return pnd.f0 + pnd.a*pnd.b*phi + pnd.b*nu

def phi_rns(phi, nu, sigma_n, f, pnd):

    F=pnd.k*(1 - spsi*np.exp(phi))*(spsi + cpsi*f_rns(phi, nu, pnd)) + sigma_n*(np.exp(phi)-np.exp(-nu))
    F=F/(pnd.a*sigma_n + pnd.eta*np.exp(phi))

    return F

def th_rns(phi,nu):

    F=np.exp(-nu)-np.exp(phi)

    return F

def sigma_rns(phi, pnd):

    F=pnd.b*pnd.k*(spsi*np.exp(phi) - 1)*cpsi

    return F

def rkf(phi,nu, sigma_n, f, h, pnd, pc):
    

    c21=1/4
    c31=3/32
    c32=9/32
    c41=1932/2197
    c42=-7200/2197
    c43=7296/2197
    c51=439/216
    c52=-8
    c53=3680/513
    c54=-845/4104
    c61=-8/27
    c62=2
    c63=-3544/2565
    c64=1859/4104
    c65=-11/40

    r1  =  1/360
    r3  = -128/4275
    r4  = -2197/75240
    r5  =  1/50
    r6  =  2/55

    c1  =  25/216
    c3  =  1408/2565
    c4  =  2197/4104
    c5  = -1/5
    c6 = 2/55

    passtep=0
    iterrk=0
    
    

    while passtep==0 and iterrk <= pc.nitrkmax:
        iterrk+=1


        # Compute values needed to compute truncation error estimate and
        # the 4th order RK estimate.
            
        
        k1 = h * phi_rns(phi,nu,sigma_n,f,pnd)
        l1 = h * th_rns(phi,nu)
        m1 = h * sigma_rns(phi,pnd)
        n1 = h * f_rns(phi,nu,pnd)


        dphi = c21*k1
        dnu = c21*l1
        dsigma = c21*m1
        df = c21*n1
        
        k2 = h * phi_rns(phi+dphi,nu+dnu,sigma_n+dsigma,f+df,pnd)
        l2 = h * th_rns(phi+dphi,nu+dnu)
        m2 = h * sigma_rns(phi+dphi,pnd)
        n2 = h * f_rns(phi+dphi,nu+dnu,pnd)

        dphi = c31*k1 + c32*k2
        dnu = c31*l1 + c32*l2
        dsigma = c31*m1 + c32*m2
        df = c31*n1 + c32*n2

        k3 = h * phi_rns(phi+dphi,nu+dnu,sigma_n+dsigma,f+df,pnd)
        l3 = h * th_rns(phi+dphi,nu+dnu)
        m3 = h * sigma_rns(phi+dphi,pnd)
        n3 = h * f_rns(phi+dphi,nu+dnu,pnd)

        dphi = c41*k1 + c42*k2 + c43*k3
        dnu = c41*l1 + c42*l2 + c43*l3
        dsigma = c41*m1 + c42*m2 + c43*m3
        df = c41*n1 + c42*n2 + c43*n3

        k4 = h * phi_rns(phi+dphi,nu+dnu,sigma_n+dsigma,f+df,pnd)
        l4 = h * th_rns(phi+dphi,nu+dnu)
        m4 = h * sigma_rns(phi+dphi,pnd)
        n4 = h * f_rns(phi+dphi,nu+dnu,pnd)

        dphi = c51*k1 + c52*k2 + c53*k3 + c54*k4
        dnu = c51*l1 + c52*l2 + c53*l3 + c54*l4
        dsigma = c51*m1 + c52*m2 + c53*m3 + c54*m4
        df = c51*n1 + c52*n2 + c53*n3 + c54*n4

        k5 = h * phi_rns(phi+dphi,nu+dnu,sigma_n+dsigma,f+df,pnd)
        l5 = h * th_rns(phi+dphi,nu+dnu)
        m5 = h * sigma_rns(phi+dphi,pnd)
        n5 = h * f_rns(phi+dphi,nu+dnu,pnd)

        dphi = c61*k1 + c62*k2 + c63*k3 + c64*k4 + c65*k5
        dnu = c61*l1 + c62*l2 + c63*l3 + c64*l4 + c65*l5
        dsigma = c61*m1 + c62*m2 + c63*m3 + c64*m4 + c65*m5
        df = c61*n1 + c62*n2 + c63*n3 + c64*n4 + c65*n5

        k6 = h * phi_rns(phi+dphi,nu+dnu,sigma_n+dsigma,f+df,pnd)
        l6 = h * th_rns(phi+dphi,nu+dnu)
        m6 = h * sigma_rns(phi+dphi,pnd)
        n6 = h * f_rns(phi+dphi,nu+dnu,pnd)


        
        # Error estimation
        err=r1*np.array([k1, l1, m1, n1])+r3*np.array([k3, l3, m3, n3])+r4*np.array([k4, l4, m4, n4])+r5*np.array([k5, l5, m5, n5])+r6*np.array([k6, l6, m6, n6])
        errmax=np.max(abs(err))

        if errmax <= pc.tol and errmax>=0:
            passtep=1
            
            phi=phi+c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 +c6 * k6
            nu=nu+c1 * l1 + c3 * l3 + c4 * l4 + c5 * l5 + c6*l6
            sigma_n=sigma_n+c1 * m1 + c3 * m3 + c4 * m4 + c5 * m5 + c6*m6
            f=f+c1 * n1 + c3 * n3 + c4 * n4 + c5 * n5 + c6*n6
            
            
            if errmax>0:
                h=np.min(np.array([pc.safe*h*((pc.tol/errmax)**0.25), pc.hmax]))


        if errmax>pc.tol:
            h=np.max(np.array([pc.safe*h*((pc.tol/errmax)**0.25), pc.hmin]))

   

    return phi, nu, sigma_n, f, dphi, dnu, dsigma, df, h # type: ignore

#-------------------------------------------#
# Iterations
#-------------------------------------------#

if __name__ == "__main__": # to allow import without running the simulation
    T=np.array([t])
    Phi=np.array([phi])
    Nu=np.array([nu])
    Sigma_n=np.array([sigma_n])
    F = np.array([f])
    Tau = np.array([f*sigma_n])
    Dphi=np.array([])
    Dnu=np.array([])

    for iter in range(0,pc.nitmax,1):
        #--update phi, nu and h
        phi, nu, sigma_n, f, dphi, dnu, dsigma, df, h = rkf(phi, nu, sigma_n, f, h, pnd, pc)

        #--update time
        t+=h

        #--store results
        T=np.append(T,[t])
        Phi=np.append(Phi,[phi])
        Nu=np.append(Nu,[nu])
        Sigma_n=np.append(Sigma_n,[sigma_n])
        F=np.append(F,[f])
        Tau = np.append(Tau, [f*sigma_n])
        Dphi=np.append(Dphi,[dphi])
        Dnu=np.append(Dnu,[dnu])


    V=np.exp(Phi)
    Dt=np.diff(T)
    Phipoint=Dphi/Dt
    Vpoint=V[1:]*Phipoint


    # save results
    r = Result(T, V, Vpoint, Nu, Phi, Phipoint, Tau, Sigma_n, pd, pnd, pc) # add filename if needed (filename = "custom_name.pkl")
    r.save_results('PR_QDYN_RNS_modele_oriente')