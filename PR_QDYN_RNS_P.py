import numpy as np
import matplotlib.pyplot as plt
import math
import pickle
from PR_QDYN_RNS import Result, ParamMec, NdParamMec, ParamComp



# -------------------------------------#
# Dimensional Mechanical parameter definition
# -------------------------------------#

pd = ParamMec(k_rigidity=3.0E10,  # rigidity (Pa)
              a_fric=0.005,  # direct effect coefficient
              b_fric=0.01,  # evolution effect coefficient
              eta_visc=1.0E18,  # viscosity (Pa.s)
              sigma_n=50.0E6,  # normal stress (Pa)
              dc= 1e-4,  # critical slip distance (m)
              V_p=1.0e-9, # tectonic speed (m/s)
              f0 = 0.6,
              r = 1.0e1,
              c = 6.8e-2,
              P0 = 25e5,
              mu = 20e9)

# -------------------------------------#
# ND Mechanical parameter definition
# -------------------------------------#
pnd = NdParamMec(a=1.3, eta=1.0E-11, k=0.4)

# pnd=NdParamMec(a = pd.a_fric/pd.b_fric, k = pd.k_rigidity*pd.dc/(pd.sigma_n*pd.b_fric), eta = pd.eta_visc*pd.V_p/(pd.b_fric*pd.sigma_n))

# -------------------------------------------#
# Computational parameter definition
# -------------------------------------------#
pc = ParamComp(tol=1.0E-10,
               nitrkmax=30,
               nitmax=10000,
               hmin=1.0E-12,
               hmax=1.0E10,
               safe=0.8)

# -------------------------------------------#
# Initial conditions (ND variables)
# -------------------------------------------#
v = 0.1  # initial normalized slip rate
th = 1 / v  # initial normalized state variable
t = 0.0001  # initial time
h = 0.001  # initial time step

phi = np.log(v)
nu = np.log(th)



def frns(phi, nu, t, pd, pnd):
    P0barre = (pd.P0/pd.sigma_n)
    rbarre = (pd.r*pd.b_fric*pd.sigma_n) / (pd.mu*pd.dc)
    cbarre = pd.c * (pd.b_fric**2 * pd.sigma_n**2)/(pd.V_p * pd.mu**2 * pd.dc)

    F = -pnd.k * (np.exp(phi)-1)
    F += (np.exp(phi) - np.exp(-nu)) * (1 - P0barre*math.erfc(rbarre / 2*math.sqrt(cbarre * t)))
    F += - (pd.f0/pd.b_fric + pnd.a*phi + nu) * (P0barre/t) * (rbarre/2*math.sqrt(np.pi*cbarre*t)) * np.exp(-(rbarre/2*math.sqrt(cbarre * t))**2)
    F /= pnd.eta*np.exp(phi) + pnd.a*(1 - P0barre*math.erfc(rbarre / 2*math.sqrt(cbarre * t)))
    F *= np.exp(phi)

    return F


def grns(phi, nu):
    G = 1/np.exp(nu) - np.exp(phi)
    return G


def rkf(phi, nu, t, h, pnd, pc):
    c21 = 1 / 4
    c31 = 3 / 32
    c32 = 9 / 32
    c41 = 1932 / 2197
    c42 = -7200 / 2197
    c43 = 7296 / 2197
    c51 = 439 / 216
    c52 = -8
    c53 = 3680 / 513
    c54 = -845 / 4104
    c61 = -8 / 27
    c62 = 2
    c63 = -3544 / 2565
    c64 = 1859 / 4104
    c65 = -11 / 40

    r1 = 1 / 360
    r3 = -128 / 4275
    r4 = -2197 / 75240
    r5 = 1 / 50
    r6 = 2 / 55

    c1 = 25 / 216
    c3 = 1408 / 2565
    c4 = 2197 / 4104
    c5 = -1 / 5
    c6 = 2 / 55

    passtep = 0
    iterrk = 0

    while passtep == 0 and iterrk <= pc.nitrkmax:
        iterrk += 1

        # Compute values needed to compute truncation error estimate and
        # the 4th order RK estimate.

        k1 = h * frns(phi, nu, t, pd, pnd)
        l1 = h * grns(phi, nu)

        dphi = c21 * k1
        dnu = c21 * l1

        k2 = h * frns(phi + dphi, nu + dnu, t, pd, pnd)
        l2 = h * grns(phi + dphi, nu + dnu)

        dphi = c31 * k1 + c32 * k2
        dnu = c31 * l1 + c32 * l2

        k3 = h * frns(phi + dphi, nu + dnu, t, pd, pnd)
        l3 = h * grns(phi + dphi, nu + dnu)

        dphi = c41 * k1 + c42 * k2 + c43 * k3
        dnu = c41 * l1 + c42 * l2 + c43 * l3

        k4 = h * frns(phi + dphi, nu + dnu, t, pd, pnd)
        l4 = h * grns(phi + dphi, nu + dnu)

        dphi = c51 * k1 + c52 * k2 + c53 * k3 + c54 * k4
        dnu = c51 * l1 + c52 * l2 + c53 * l3 + c54 * l4

        k5 = h * frns(phi + dphi, nu + dnu, t, pd, pnd)
        l5 = h * grns(phi + dphi, nu + dnu)

        dphi = c61 * k1 + c62 * k2 + c63 * k3 + c64 * k4 + c65 * k5
        dnu = c61 * l1 + c62 * l2 + c63 * l3 + c64 * l4 + c65 * l5

        k6 = h * frns(phi + dphi, nu + dnu, t, pd, pnd)
        l6 = h * grns(phi + dphi, nu + dnu)

        # Error estimation
        err = r1 * np.array([k1, l1]) + r3 * np.array([k3, l3]) + r4 * np.array([k4, l4]) + r5 * np.array(
            [k5, l5]) + r6 * np.array([k6, l6])
        errmax = np.max(abs(err))

        if errmax <= pc.tol and errmax >= 0:
            passtep = 1

            phi = phi + c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6
            nu = nu + c1 * l1 + c3 * l3 + c4 * l4 + c5 * l5 + c6 * l6

            if errmax > 0:
                h = np.min(np.array([pc.safe * h * ((pc.tol / errmax) ** 0.25), pc.hmax]))

        if errmax > pc.tol:
            h = np.max(np.array([pc.safe * h * ((pc.tol / errmax) ** 0.25), pc.hmin]))

    return phi, nu, dphi, dnu, h  # type: ignore


# -------------------------------------------#
# Iterations
# -------------------------------------------#

if __name__ == "__main__":  # to allow import without running the simulation
    T = np.array([t])
    Dt = np.array([])
    Phi = np.array([phi])
    Nu = np.array([nu])
    Dphi = np.array([])
    Dnu = np.array([])

    for iter in range(0, pc.nitmax, 1):
        # --update phi, nu and h
        phi, nu, dphi, dnu, h = rkf(phi, nu, t, h, pnd, pc)

        # --update time
        t += h

        # --store results
        T = np.append(T, [t])
        Dt = np.append(Dt, [h])
        Phi = np.append(Phi, [phi])
        Nu = np.append(Nu, [nu])
        Dphi = np.append(Dphi, [dphi])
        Dnu = np.append(Dnu, [dnu])

    V = np.exp(Phi)
    Vln = np.log(V)
    Phipoint = Dphi / Dt
    Vpoint = V[1:] * Phipoint
    #Vpointln = np.log(Vpoint)

    # save results
    Result(T, V, Vpoint, Nu, Phi, pd, pnd, pc, filename="testPression.pkl").save_results()  # add filename if needed (filename = "custom_name.pkl")
