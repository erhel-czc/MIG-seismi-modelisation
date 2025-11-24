import os
import numpy as np
import matplotlib
# # Use a non-interactive backend when no DISPLAY is available (headless environments)
# if not os.environ.get('DISPLAY'):
#     matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle

class ParamMec:
    "Dimensional Mechanical parameters"

    def __init__(self, k_rigidity, a_fric, b_fric, eta_visc, sigma_n, dc, V_p):
        self.k_rigidity=k_rigidity
        self.a_fric=a_fric
        self.b_fric=b_fric
        self.eta_visc=eta_visc
        self.sigma_n=sigma_n
        self.dc=dc
        self.V_p=V_p

        

class NdParamMec:
    "Non dimensional Mechanical parameters"

    def __init__(self, a, eta, k):
        self.a=a
        self.eta=eta
        self.k=k

class ParamComp:
    "Computational parameters"

    def __init__(self, tol, nitrkmax, nitmax, hmin, hmax, safe):
        self.tol=tol # error tolerance
        self.nitrkmax=nitrkmax # maximum number of iteration in a rkf step
        self.nitmax=nitmax # maximum number of iterations
        self.hmin=hmin # minimum time step
        self.hmax=hmax # maximum time step (CFL for diffusion equation)
        self.safe=safe

class Result:
    "Results storage"

    def __init__(self, T, V, Vpoint, Nu, Phi, pd, pnd, pc, filename = ''):
        self.T=T
        self.V=V
        self.Vpoint=Vpoint
        self.Nu=Nu
        self.Phi=Phi
        self.pd=pd
        self.pnd=pnd
        self.pc=pc
        if filename=='':
            self.filename= f"{pnd.a}_{pnd.eta}_{pnd.k}.pkl"
        else:
            self.filename=filename

    def save_results(self):
        with open(f"Results/{self.filename}", 'wb') as f:
            pickle.dump(self, f)
        f.close()

    @staticmethod
    def load_results(filename):
        with open(f"Results/{filename}", 'rb') as f:
            data = pickle.load(f)
        f.close()
        return data

    def slip_rate_evolution(self):
        plt.figure('Slip rate evolution')

        plt.plot(self.T, np.log(self.V), '-k')
        plt.xlabel('Time (ND)')
        plt.ylabel('Log Slip rate (ND)')
        plt.title(r'Slip rate evolution ($\kappa$=%.2f, $\alpha$=%.2f)' % (self.pnd.k, self.pnd.a))
        plt.grid()

        # plt.xlim([0, 100])
        plt.savefig('images/modif_parameters/slip_rate_k%.2f_a%.2f.pdf' % (self.pnd.k, self.pnd.a))

    def phase_portrait(self):
        plt.figure('Phase portrait')

        sc = plt.scatter(self.V[1:], self.Vpoint, c=self.T[1:], cmap='viridis', marker='+')
        plt.plot(self.V[1:], self.Vpoint, '-k', alpha=0.3)
        plt.xlabel('Speed (ND)')
        plt.ylabel('Acceleration (ND)')
        plt.title(r'Phase portrait ($\kappa$=%.2f, $\alpha$=%.2f)' % (self.pnd.k, self.pnd.a))
        plt.colorbar(sc, label='Time progression')
        plt.grid()

        plt.savefig('images/modif_parameters/phase_portrait_k%.2f_a%.2f.pdf' % (self.pnd.k, self.pnd.a))

        plt.show()




#-------------------------------------#
# Dimensional Mechanical parameter definition
#-------------------------------------#

pd = ParamMec(k_rigidity=3.0E10, # rigidity (Pa)
              a_fric=0.005,     # direct effect coefficient
              b_fric=0.01,      # evolution effect coefficient
              eta_visc=1.0E18,  # viscosity (Pa.s)
              sigma_n=50.0E6,   # normal stress (Pa)
              dc=0.01,          # critical slip distance (m)
              V_p=1.0E-6)       # tectonic speed (m/s)

#-------------------------------------#
# ND Mechanical parameter definition
#-------------------------------------#
pnd=NdParamMec(a=0.5, eta=1.0E-11, k=0.4)

#pnd=NdParamMec(a = pd.a_fric/pd.b_fric, k = pd.k_rigidity*pd.dc/(pd.sigma_n*pd.b_fric), eta = pd.eta_visc*pd.V_p/(pd.b_fric*pd.sigma_n))

#-------------------------------------------#
# Computational parameter definition
#-------------------------------------------#
pc=ParamComp(tol=1.0E-10,
             nitrkmax=30,
             nitmax=10000,
             hmin=1.0E-12,
             hmax=1.0E10,
             safe=0.8)

#-------------------------------------------#
# Initial conditions (ND variables)
#-------------------------------------------#
v=0.1       # initial normalized slip rate
th=1/v      # initial normalized state variable
t=0.0       # initial time
h=0.001     # initial time step

phi=np.log(v)
nu=np.log(th)

def frns(phi,nu,pnd):

    F=-pnd.k*(np.exp(phi)-1)+np.exp(phi)-np.exp(-nu)
    F=F/(pnd.a+pnd.eta*np.exp(phi))

    return F

def grns(phi,nu):

    F=np.exp(-nu)-np.exp(phi)

    return F

def rkf(phi,nu,h,pnd,pc):
    

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
            
        
        k1 = h * frns(phi,nu,pnd)
        l1 = h * grns(phi,nu)


        dphi = c21*k1
        dnu = c21*l1
        
        k2 = h * frns(phi+dphi,nu+dnu,pnd)
        l2 = h * grns(phi+dphi,nu+dnu)


        dphi = c31*k1+c32*k2
        dnu = c31*l1+c32*l2
        
        k3 = h * frns(phi+dphi,nu+dnu,pnd)
        l3 = h * grns(phi+dphi,nu+dnu)


        dphi = c41*k1+c42*k2+c43*k3
        dnu = c41*l1+c42*l2+c43*l3
        
        k4 = h * frns(phi+dphi,nu+dnu,pnd)
        l4 = h * grns(phi+dphi,nu+dnu)



        dphi = c51*k1+c52*k2+c53*k3+c54*k4
        dnu = c51*l1+c52*l2+c53*l3+c54*l4

        k5 = h * frns(phi+dphi,nu+dnu,pnd)
        l5 = h * grns(phi+dphi,nu+dnu)


        dphi = c61*k1+c62*k2+c63*k3+c64*k4+c65*k5
        dnu = c61*l1+c62*l2+c63*l3+c64*l4+c65*l5
        
        k6 = h * frns(phi+dphi,nu+dnu,pnd)
        l6 = h * grns(phi+dphi,nu+dnu)


        
        # Error estimation
        err=r1*np.array([k1, l1])+r3*np.array([k3, l3])+r4*np.array([k4, l4])+r5*np.array([k5, l5])+r6*np.array([k6, l6])
        errmax=np.max(abs(err))

        if errmax <= pc.tol and errmax>=0:
            passtep=1
            
            phi=phi+c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 +c6 * k6
            nu=nu+c1 * l1 + c3 * l3 + c4 * l4 + c5 * l5 + c6*l6
            
            
            if errmax>0:
                h=np.min(np.array([pc.safe*h*((pc.tol/errmax)**0.25), pc.hmax]))


        if errmax>pc.tol:
            h=np.max(np.array([pc.safe*h*((pc.tol/errmax)**0.25), pc.hmin]))

   

    return phi,nu,dphi,dnu,h # type: ignore

#-------------------------------------------#
# Iterations
#-------------------------------------------#

if __name__ == "__main__": #to allow import without running the simulation
    T=np.array([t])
    Phi=np.array([phi])
    Nu=np.array([nu])
    Dphi=np.array([])
    Dnu=np.array([])

    for iter in range(0,pc.nitmax,1):
        #--update phi, nu and h
        phi,nu,dphi,dnu,h = rkf(phi,nu,h,pnd,pc)

        #--update time
        t+=h

        #--store results
        T=np.append(T,[t])
        Phi=np.append(Phi,[phi])
        Nu=np.append(Nu,[nu])
        Dphi=np.append(Dphi,[dphi])
        Dnu=np.append(Dnu,[dnu])


    V=np.exp(Phi)
    Vln=np.log(V)
    Dt=np.diff(T)
    Phipoint=Dphi/Dt
    Vpoint=V[1:]*Phipoint
    Vpointln=np.log(Vpoint)

    #save results
    Result(T, V, Vpoint, Nu, Phi, pd, pnd, pc).save_results() #add filename if needed (filename = "custom_name.pkl")



#-------------------------------------------#
# Plot
#-------------------------------------------#
#
# """ Slip rate evolution """
# plt.figure('Slip rate evolution')
#
# plt.plot(T,Vln,'-k')
# plt.xlabel('Time (ND)')
# plt.ylabel('Log Slip rate (ND)')
# plt.title(r'Slip rate evolution ($\kappa$=%.2f, $\alpha$=%.2f)' % (pnd.k, pnd.a))
# plt.grid()
#
# #plt.xlim([0, 100])
# plt.savefig('images/modif_parameters/slip_rate_k%.2f_a%.2f.pdf' % (pnd.k, pnd.a))
#
# """ Phase portrait """
# plt.figure('Phase portrait')
#
# sc = plt.scatter(V[1:], Vpoint, c=T[1:], cmap='viridis', marker='+')
# plt.plot(V[1:], Vpoint, '-k', alpha=0.3)
# plt.xlabel('Speed (ND)')
# plt.ylabel('Acceleration (ND)')
# plt.title(r'Phase portrait ($\kappa$=%.2f, $\alpha$=%.2f)' % (pnd.k, pnd.a))
# plt.colorbar(sc, label='Time progression')
# plt.grid()
#
# plt.savefig('images/modif_parameters/phase_portrait_k%.2f_a%.2f.pdf' % (pnd.k, pnd.a))
#
# plt.show()
#
# #plt.plot(T,np.log10(np.exp(Nu)),'-+k')
# #plt.show()
