import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom
from pathlib import Path


def calcShearCart(t, x, pars):
    rho0 = 1.0
    gam = pars['Adiabatic_Index']
    nu = pars['Viscosity']
    x0 = pars['Init_Par1']
    cs0 = pars['Init_Par2']
    sig0 = pars['Init_Par3']
    v0 = pars['Init_Par4']

    rho = np.empty(x.shape)
    rho[:] = rho0

    if pars['Isothermal']:
        cs2 = 1.0 / pars['Mach_Number']**2
        P = rho * cs2 / gam
    else:
        P = rho * cs0*cs0 / gam

    vx = np.zeros(x.shape)

    t0 = sig0*sig0 / (2*nu)
    sig = np.sqrt(2*nu*(t+t0))

    vy = v0 * np.exp(-0.5*(x-x0)*(x-x0) / (sig*sig)) * sig0/sig

    return rho, P, vx, vy


def analyzeSingle(filename):

    opts = util.loadOpts(filename)
    pars = util.loadPars(filename)
    print("Loading " + filename)
    t, x1, x2, x3, prim, dat = util.loadCheckpoint(filename)

    nx = pars['Num_R']
    v0 = pars['Init_Par4']

    rho = prim[:, 0]
    P = prim[:, 1]
    v1 = prim[:, 2]
    v2 = prim[:, 3]
    v3 = prim[:, 4]

    x, y, z = geom.getXYZ(x1, x2, x3, opts, pars)
    vx, vy, vz = geom.getVXYZ(x1, x2, x3, v1, v2, v3, opts)

    rhoS, PS, vxS, vyS = calcShearCart(t, x, pars)

    dV = geom.getDV(dat, opts, pars)

    maskDistance = 0.0

    dVmask = dV.copy()
    # dVmask[(x1 < pars['R_Min']+maskDistance)] = 0.0
    dVmask[(x1 > pars['R_Max']-maskDistance)] = 0.0

    errRho = geom.integrate(np.fabs(rho-rhoS), dat, opts, pars, dVmask)
    errP = geom.integrate(np.fabs(P-PS), dat, opts, pars, dVmask)
    errVx = geom.integrate(np.fabs(vx-vxS), dat, opts, pars, dVmask)
    errVy = geom.integrate(np.fabs(vy-vyS), dat, opts, pars, dVmask)

    X = np.linspace(pars['R_Min'], pars['R_Max'], 500)
    _, _, _, VY = calcShearCart(t, X, pars)

    name = (Path(filename).stem).split('_')[-1]
    figname = 'shearCart_vy_{0:s}.png'.format(name)

    fig, ax = plt.subplots(1, 1)
    ax.plot(x, vy, ls='', marker='+', color='k')
    ax.plot(X, VY)
    ax.set_xlim(pars['R_Min'], pars['R_Max'])
    ax.set_ylim(-0.2 * v0, 1.2*v0)
    print("Saving " + figname)
    fig.savefig(figname)
    plt.close(fig)

    return t, nx, errRho, errP, errVx, errVy


def analyze(filenames):

    N = len(filenames)

    t = np.empty(N)
    nx = np.empty(N)
    errRho = np.empty(N)
    errVx = np.empty(N)
    errVy = np.empty(N)
    errVz = np.empty(N)
    errP = np.empty(N)
    errS = np.empty(N)
    errJp = np.empty(N)
    errJm = np.empty(N)

    for i, f in enumerate(filenames):
        dat = analyzeSingle(f)
        t[i] = dat[0]
        nx[i] = dat[1]
        errRho[i] = dat[2]
        errP[i] = dat[3]
        errVx[i] = dat[4]
        errVy[i] = dat[5]

    makeErrPlot(t, nx, errRho, "acousticwave_errRho", r"$L_1(\rho)$")
    makeErrPlot(t, nx, errP, "acousticwave_errP", r"$L_1(P)$")
    makeErrPlot(t, nx, errVx, "acousticwave_errVx", r"$L_1(v_x)$")
    makeErrPlot(t, nx, errVy, "acousticwave_errVy", r"$L_1(v_y)$")


def makeErrPlot(t, nx, err, name, label):
    print(name)
    print(err)
    T = np.unique(t)
    NX = np.unique(nx)
    figname = name + "_nx.png"
    if len(NX) > 1:
        fig, ax = plt.subplots(1, 1)
        for tt in T:
            ind = t == tt
            ax.plot(nx[ind], err[ind], marker='+', ms=10, mew=1, ls='')
        nn = np.logspace(math.log10(NX[1]), math.log10(NX[-1]),
                         num=10, base=10.0)
        ax.plot(nn, err.max() * np.power(nn/nn[0], -2), ls='--', lw=2,
                color='grey')
        ax.set_xlabel(r"$n_x$")
        ax.set_ylabel(label)
        ax.set_xscale('log')
        if (err > 0).any():
            ax.set_yscale('log')
            ylim = ax.get_ylim()
            if ylim[1] < 10*ylim[0]:
                ax.set_ylim(0.25*ylim[0], 4*ylim[1])
        print("Saving " + figname)
        fig.savefig(figname)
        plt.close(fig)

    figname = name + "_t.png"
    if len(T) > 1:
        fig, ax = plt.subplots(1, 1)
        for nn in NX:
            ind = nn == nx
            ax.plot(t[ind], err[ind])
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(label)
        ax.set_xscale('log')
        if (err > 0).any():
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
