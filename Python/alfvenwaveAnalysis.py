import sys
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discoUtil as du
import discoGeom as dg

def psiProf(x, pars):

    choice = int(pars['Init_Par0'])

    if choice == 0:
        lam = pars['Init_Par4']
        phase0 = pars['Init_Par7']
        psi = np.pi*phase0 + 2*np.pi*x/lam
    elif choice == 1:
        lam = pars['Init_Par4']
        x0 = pars['Init_Par7']
        L = pars['Init_Par8']
        psi = 2*np.pi*L/lam * 0.5*(1.0 + np.tanh((x-x0)/L))
    else:
        psi = 0.0

    return psi

def calcAlfvenWave(t, x, pars):
    rho0 = 1.0
    cs0 = 1.0
    gam = pars['Adiabatic_Index']
    v0 = pars['Init_Par1']
    cA = pars['Init_Par2']
    a = pars['Init_Par3']
    lam = pars['Init_Par4']

    xi = x - cA*t

    psi = psiProf(xi, pars)
    cp = np.cos(psi)
    sp = np.sin(psi)
    B0 = math.sqrt(rho0 * cA*cA)

    rho = np.empty(x.shape)
    P = np.empty(x.shape)
    vx = np.empty(x.shape)
    vy = np.empty(x.shape)
    vz = np.empty(x.shape)
    Bx = np.empty(x.shape)
    By = np.empty(x.shape)
    Bz = np.empty(x.shape)

    rho[:] = rho0
    P[:] = rho0*cs0*cs0/gam
    vx[:] = v0
    vy[:] = -B0*a*cp/math.sqrt(rho0)
    vz[:] = -B0*a*sp/math.sqrt(rho0)
    Bx[:] = B0
    By[:] = B0*a*cp
    Bz[:] = B0*a*sp

    return rho, P, vx, vy, vz, Bx, By, Bz


def analyzeSingle(filename):

    opts = du.loadOpts(filename)
    pars = du.loadPars(filename)
    print("Loading " + filename)
    t, x1, x2, x3, prim, dat = du.loadCheckpoint(filename)

    nx = pars['Num_R']
    gam = pars['Adiabatic_Index']

    rho = prim[:,0]
    P = prim[:,1]
    v1 = prim[:,2]
    v2 = prim[:,3]
    v3 = prim[:,4]
    B1 = prim[:,5]
    B2 = prim[:,6]
    B3 = prim[:,7]
    eS = P * np.power(rho, -gam)

    x, y, z = dg.getXYZ(x1, x2, x3, opts, pars)
    vx, vy, vz = dg.getVXYZ(x1, x2, x3, v1, v2, v3, opts)
    Bx, By, Bz = dg.getBXYZ(x1, x2, x3, B1, B2, B3, opts)

    gam = pars['Adiabatic_Index']
    phi0 = pars['Init_Par5'] * np.pi
    ct = pars['Init_Par6']
    st = math.sqrt((1.0 - ct)*(1.0+ct))
    cp = math.cos(phi0)
    sp = math.sin(phi0)
    
    ix, iy, iz = st*cp, st*sp, ct
    jx, jy, jz = -sp, cp, 0.0
    kx, ky, kz = -ct*cp, -ct*sp, st

    X = ix*x + iy*y + iz*z
    #VX = ix*vx + iy*vy + iz*vz
    #VY = jx*vx + jy*vy + jz*vz
    #VZ = kx*vx + ky*vy + kz*vz
    #BX = ix*Bx + iy*By + iz*Bz
    #BY = jx*Bx + jy*By + jz*Bz
    #BZ = kx*Bx + ky*By + kz*Bz

    rhoS, PS, ViS, VjS, VkS, BiS, BjS, BkS = calcAlfvenWave(t, X, pars)
    eSS = PS * np.power(rhoS, -gam)

    dV = dg.getDV(dat, opts, pars)

    errRho = dg.integrate(np.fabs(rho-rhoS), dat, opts, pars, dV)
    errP = dg.integrate(np.fabs(P-PS), dat, opts, pars, dV)
    errVx = dg.integrate(np.fabs(vx-(ix*ViS+jx*VjS+kx*VkS)), dat,opts,pars,dV)
    errVy = dg.integrate(np.fabs(vy-(iy*ViS+jy*VjS+ky*VkS)), dat,opts,pars,dV)
    errVz = dg.integrate(np.fabs(vz-(iz*ViS+jz*VjS+kz*VkS)), dat,opts,pars,dV)
    errBx = dg.integrate(np.fabs(Bx-(ix*BiS+jx*BjS+kx*BkS)), dat,opts,pars,dV)
    errBy = dg.integrate(np.fabs(By-(iy*BiS+jy*BjS+ky*BkS)), dat,opts,pars,dV)
    errBz = dg.integrate(np.fabs(Bz-(iz*BiS+jz*BjS+kz*BkS)), dat,opts,pars,dV)
    errS = dg.integrate(np.fabs(eS-eSS), dat, opts, pars, dV)

    return t, nx, errRho, errP, errVx, errVy, errVz, errBx, errBy, errBz, errS

def analyze(filenames):

    N = len(filenames)

    t = np.empty(N)
    nx = np.empty(N)
    errRho = np.empty(N)
    errP = np.empty(N)
    errVx = np.empty(N)
    errVy = np.empty(N)
    errVz = np.empty(N)
    errBx = np.empty(N)
    errBy = np.empty(N)
    errBz = np.empty(N)
    errS = np.empty(N)

    for i,f in enumerate(filenames):
        dat = analyzeSingle(f)
        t[i] = dat[0]
        nx[i] = dat[1]
        errRho[i] = dat[2]
        errP[i] = dat[3]
        errVx[i] = dat[4]
        errVy[i] = dat[5]
        errVz[i] = dat[6]
        errBx[i] = dat[7]
        errBy[i] = dat[8]
        errBz[i] = dat[9]
        errS[i] = dat[10]

    makeErrPlot(t, nx, errRho, "alfvenwave_errRho", r"$L_1(\rho)$")
    makeErrPlot(t, nx, errP, "alfvenwave_errP", r"$L_1(P)$")
    makeErrPlot(t, nx, errVx, "alfvenwave_errVx", r"$L_1(v_x)$")
    makeErrPlot(t, nx, errVy, "alfvenwave_errVy", r"$L_1(v_y)$")
    makeErrPlot(t, nx, errVz, "alfvenwave_errVz", r"$L_1(v_z)$")
    makeErrPlot(t, nx, errBx, "alfvenwave_errBx", r"$L_1(B_x)$")
    makeErrPlot(t, nx, errBy, "alfvenwave_errBy", r"$L_1(B_y)$")
    makeErrPlot(t, nx, errBz, "alfvenwave_errBz", r"$L_1(B_z)$")
    makeErrPlot(t, nx, errS, "alfvenwave_errS", r"$L_1(P/\rho^\gamma)$")


def makeErrPlot(t, nx, err, name, label):
    print(name)
    print(err)
    T = np.unique(t)
    NX = np.unique(nx)
    figname = name + "_nx.png"
    if len(NX) > 1:
        fig, ax = plt.subplots(1,1)
        for tt in T:
            ind = t==tt
            ax.plot(nx[ind], err[ind], marker='+', ms=10, mew=1, ls='')
        nn = np.logspace(math.log10(NX[1]), math.log10(NX[-1]), 
                            num=10, base=10.0)
        ax.plot(nn, err.max() * np.power(nn/nn[0], -2), ls='--', lw=2, 
                    color='grey')
        ax.set_xlabel(r"$n_x$")
        ax.set_ylabel(label)
        ax.set_xscale('log')
        if (err>0).any():
            ax.set_yscale('log')
            ylim = ax.get_ylim()
            if ylim[1] < 10*ylim[0]:
                ax.set_ylim(0.25*ylim[0], 4*ylim[1])
        print("Saving " + figname)
        fig.savefig(figname)
        plt.close(fig)

    figname = name + "_t.png"
    if len(T) > 1:
        fig, ax = plt.subplots(1,1)
        for nn in NX:
            ind = nn==nx
            ax.plot(t[ind], err[ind])
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(label)
        ax.set_xscale('log')
        if (err>0).any():
            ax.set_yscale('log')
            ylim = ax.get_ylim()
            if ylim[1] < 10*ylim[0]:
                ax.set_ylim(0.25*ylim[0], 4*ylim[1])
        print("Saving " + figname)
        fig.savefig(figname)
        plt.close(fig)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need a checkpoint dude!")
        sys.exit()

    filenames = sys.argv[1:]
    analyze(filenames)
