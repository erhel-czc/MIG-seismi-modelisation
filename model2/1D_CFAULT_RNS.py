import numpy as np
import math
import matplotlib.pyplot as plt
import numpy.random
import result1D as r

#-------------------------------------------#
# Computational parameter definition
#-------------------------------------------#
dx=0.2                           # grid size
n=128                             # number of points
tol=1.0E-7                        # error tolerance
######################
# préférer tol=1.0E-12 pour les simulations
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

class NdParamMec:
    "Non dimensional Mechanical parameters"

    def __init__(self, alpha, eta, a, b, taupinf, h):
        self.alpha = alpha
        self.eta = eta
        self.a = a
        self.b = b
        self.taupinf = taupinf
        self.h = h

class ParamComp:
    "Computational parameters"

    def __init__(self, dx, n, tol, nitrkmax, nitmax, hmin, hmax, safe):
        self.dx = dx
        self.n = n
        self.tol = tol
        self.nitrkmax = nitrkmax
        self.nitmax = nitmax
        self.hmin = hmin
        self.hmax = hmax
        self.safe = safe


#-------------------------------------------#
# Computational parameter definition
#-------------------------------------------#
pc=ParamComp(dx, n, tol, nitrkmax, nitmax, hmin, hmax, safe)

#-------------------------------------#
# ND Mechanical parameter definition
#-------------------------------------#
pnd=NdParamMec(alpha, eta,
               a=np.ones(pc.n),
               b=np.ones(pc.n),
               taupinf=0.0*np.ones(pc.n),
               h=10*pc.n*pc.dx)

x=np.arange(-0.5*pc.n*pc.dx,0.5*pc.n*pc.dx,pc.dx)
iasp=np.where(np.abs(x)<5)
ianti=np.where(np.abs(x)>=5)


phi=np.log(v0)*np.ones(pc.n)
nu=np.log(th0)*np.ones(pc.n)

phi[iasp]=np.log(0.001)
#pnd.a[ianti]=7


def gthilb(phi,pc):
    
    F=np.fft.fft(np.exp(phi))
    
  #  print(np.shape(F))
    
    k=np.fft.fftfreq(pc.n,pc.dx)
    
  #  print(np.shape(k))
    inz=np.where(np.abs(k)>0)
    iz=np.where(k==0)
   
    F[iz]=np.zeros(np.shape(iz))

    F[inz]=np.abs(k[inz]*2*math.pi)*((1+np.exp(-4*pnd.h*np.abs(k[inz]*2*math.pi)))/(1-np.exp(-4*pnd.h*np.abs(k[inz]*2*math.pi))))*F[inz]
    
    #print(np.shape(F))

    gthc=np.fft.ifft(F)
    
    gth=gthc.real
    
   # print(np.shape(gth))
    
    #gth=np.reshape(gth,(pc.n,1))
    
    return gth

def frns(phi,nu,pnd,pc):
    
    gth=gthilb(phi, pc)
    
    #print(np.shape(gth))
    #print(np.shape(pnd.taupinf))
    #print(np.shape(pnd.b))
    #print(np.shape(phi))
    #print(np.shape(nu))

    F=-0.5*gth -(np.exp(phi)-1)/pnd.h + pnd.taupinf +pnd.b*(np.exp(phi)-np.exp(-nu))
    
    #print(np.shape(F))
 
    F=F/(pnd.alpha*pnd.a+pnd.eta*np.exp(phi))
    
    #print(np.shape(F))  
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
            
        
        k1 = h * frns(phi,nu,pnd,pc)
        l1 = h * grns(phi,nu)
        
        #print(np.shape(k1))
        #print(np.shape(l1))
    


        dphi = c21*k1
        dnu = c21*l1
        
        k2 = h * frns(phi+dphi,nu+dnu,pnd,pc)
        l2 = h * grns(phi+dphi,nu+dnu)


        dphi = c31*k1+c32*k2
        dnu = c31*l1+c32*l2
        
        k3 = h * frns(phi+dphi,nu+dnu,pnd,pc)
        l3 = h * grns(phi+dphi,nu+dnu)


        dphi = c41*k1+c42*k2+c43*k3
        dnu = c41*l1+c42*l2+c43*l3
        
        k4 = h * frns(phi+dphi,nu+dnu,pnd,pc)
        l4 = h * grns(phi+dphi,nu+dnu)



        dphi = c51*k1+c52*k2+c53*k3+c54*k4
        dnu = c51*l1+c52*l2+c53*l3+c54*l4

        k5 = h * frns(phi+dphi,nu+dnu,pnd,pc)
        l5 = h * grns(phi+dphi,nu+dnu)


        dphi = c61*k1+c62*k2+c63*k3+c64*k4+c65*k5
        dnu = c61*l1+c62*l2+c63*l3+c64*l4+c65*l5
        
        k6 = h * frns(phi+dphi,nu+dnu,pnd,pc)
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

   

    return phi,nu,h

#-------------------------------------------#
# Iterations
#-------------------------------------------#
T=np.array([t])
Phi=np.array([phi])
Nu=np.array([nu])
phimax=np.array([np.max(phi)])

phiref=np.max(phi)

for iter in range(0,pc.nitmax,1):
    
    #--counter
    if iter%100 == 0:
        print("iteration ",iter," max iteration: ",pc.nitmax," time (nd): ",t," slip rate (nd): ",np.exp(np.max(phi)))
    
    #--update phi, nu and h
    phi,nu,h = rkf(phi,nu,h,pnd,pc)
    
    #--update time
    t+=h

    #--plot & store results
    if np.abs(np.max(phi)-phiref)>1:
        #--plot profiles
        phiref=np.max(phi)
        #plt.plot(x,phi)
        
        #--store
        T=np.append(T,[t])
        Phi=np.append(Phi,phi)
        Nu=np.append(Nu,nu)
        phimax=np.append(phimax,[np.max(phi)])
        
        Phi=np.reshape(Phi,(int(len(Phi)/(pc.n)),pc.n))
        Nu=np.reshape(Nu,(int(len(Nu)/(pc.n)),pc.n))


sim = r.Result1D(T=T, V=np.exp(Phi), Nu=Nu, Phi=Phi, pd=pnd, pnd=pnd, pc=pc, filename='test_simulation.pkl')
sim.save_results(folder_name='test')

#-------------------------------------------#
# Plot
#-------------------------------------------#
"""plt.contour(x,np.log10(np.max(T)-T[0:1000]),np.log10(np.exp(Phi[0:1000])),30)
plt.contour(x,T,np.log10(np.exp(Phi)),30)
plt.contour(x,np.log10(T[1:pc.nitmax+1]),np.log10(np.exp(Phi[1:pc.nitmax+1])),30)
plt.colorbar(orientation='vertical')
plt.show()

plt.plot(np.log10(np.max(T)-T[0:1000]),np.log10(np.exp(phimax[0:1000])),'-+k')
plt.xlim([-5,-4.8])
plt.plot(np.log10(T[1:pc.nitmax+1]),np.log10(np.exp(phimax[1:pc.nitmax+1])),'-+k')
plt.plot(T,np.log10(np.exp(phimax)),'-+k')
plt.show()

# plt.figure('Phi evolution')
# plt.plot(T,pnd.a*Phi+Nu,'-+k')
plt.show()

"""