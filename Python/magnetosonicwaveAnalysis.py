import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom
import calc as ca


def rhoProf(x, a, x0, L, rho0):

    rho = np.empty(x.shape)
    rho[:] = rho0
    inn = np.fabs(x-x0) < L
    xin = x[inn]
    f = 1.0 - (xin-x0)*(xin-x0)/(L*L)
    rho[inn] *= 1.0 + a*f*f*f*f

    return rho


def calcMagnetosonicWave(t, x, pars, tol=1.0e-6):
    rho0 = 1.0
    cs0 = 1.0
    gam = pars['Adiabatic_Index']
    v0 = pars['Init_Par1']
    cA0 = pars['Init_Par2']
    a = pars['Init_Par3']
    L = pars['Init_Par4']
    x0 = pars['Init_Par7']

    xR = x.max()
    xL = x.min()

    rhomax = rho0 * (1+a)
    csmax = cs0*np.power(rhomax/rho0, 0.5*(gam-1))
    cAmax = cA0*np.sqrt(rhomax/rho0)
    vmax = v0 + ca.magnetosonic_cf_int_newt(np.array([rhomax]), rho0, cs0,
                                            cA0, gam)[0]
    cfmax = math.sqrt(csmax*csmax + cAmax*cAmax)
    cf0 = math.sqrt(cs0*cs0 + cA0*cA0)

    xa = x - (1.1*(vmax+cfmax))*t
    xb = x - (0.9*(v0+cf0))*t

    err = np.empty(x.shape)
    err[:] = np.inf

    while np.fabs(err).min() > tol:
        xc = 0.5*(xa+xb)
        rho = rhoProf(xc, a, x0, L, rho0)
        cs2 = cs0*cs0*np.power(rho/rho0, gam-1)
        cA2 = cA0*cA0*rho/rho0

        cf = np.sqrt(cs2+cA2)
        v = v0 + ca.magnetosonic_cf_int_newt(rho, rho0, cs0, cA0, gam)
        xf = xc + (v+cf)*t

        err = (xf-x) / (xR-xL)
        over = err > 0
        under = err < 0
        xb[over] = xc[over]
        xa[under] = xc[under]

    xc = 0.5*(xa+xb)
    rho = rhoProf(xc, a, x0, L, rho0)
    P = rho0*cs0*cs0/gam * np.power(rho/rho0, gam)
    v = v0 + ca.magnetosonic_cf_int_newt(rho, rho0, cs0, cA0, gam)
    cA = cA0 * np.sqrt(rho/rho0)
    B = np.sqrt(rho) * cA

    return rho, P, v, B


def analyzeSingle(filename):

    opts = util.loadOpts(filename)
    pars = util.loadPars(filename)
    print("Loading " + filename)
    t, x1, x2, x3, prim, dat = util.loadCheckpoint(filename)

    nx = pars['Num_R']
    gam = pars['Adiabatic_Index']

    rho = prim[:, 0]
    P = prim[:, 1]
    v1 = prim[:, 2]
    v2 = prim[:, 3]
    v3 = prim[:, 4]
    B1 = prim[:, 5]
    B2 = prim[:, 6]
    B3 = prim[:, 7]
    eS = P * np.power(rho, -gam)

    x, y, z = geom.getXYZ(x1, x2, x3, opts, pars)
    vx, vy, vz = geom.getVXYZ(x1, x2, x3, v1, v2, v3, opts)
    Bx, By, Bz = geom.getBXYZ(x1, x2, x3, B1, B2, B3, opts)

    gam = pars['Adiabatic_Index']
    phi0 = pars['Init_Par5'] * np.pi
    ct = pars['Init_Par6']
    st = math.sqrt((1.0 - ct)*(1.0+ct))
    cp = math.cos(phi0)
    sp = math.sin(phi0)
    phiB = pars['Init_Par8'] * math.pi

    ix, iy, iz = st*cp, st*sp, ct
    jx, jy, jz = -sp, cp, 0.0
    kx, ky, kz = -ct*cp, -ct*sp, st

    X = ix*x + iy*y + iz*z

    # Jm = V + 2*cs/(gam-1)
    # Jp = V - 2*cs/(gam-1)

    rhoS, PS, vS, BS = calcMagnetosonicWave(t, X, pars, 1.0e-13)
    eSS = PS * np.power(rhoS, -gam)
    BjS = math.cos(phiB) * BS
    BkS = math.sin(phiB) * BS
    # JmS = VS + 2*csS/(gam-1)
    # JpS = VS - 2*csS/(gam-1)

    dV = geom.getDV(dat, opts, pars)

    errRho = geom.integrate(np.fabs(rho-rhoS), dat, opts, pars, dV)
    errP = geom.integrate(np.fabs(P-PS), dat, opts, pars, dV)
    errVx = geom.integrate(np.fabs(vx-ix*vS), dat, opts, pars, dV)
    errVy = geom.integrate(np.fabs(vy-iy*vS), dat, opts, pars, dV)
    errVz = geom.integrate(np.fabs(vz-iz*vS), dat, opts, pars, dV)
    errBx = geom.integrate(np.fabs(Bx-(jx*BjS+kx*BkS)), dat, opts, pars, dV)
    errBy = geom.integrate(np.fabs(By-(jy*BjS+ky*BkS)), dat, opts, pars, dV)
    errBz = geom.integrate(np.fabs(Bz-(jz*BjS+kz*BkS)), dat, opts, pars, dV)
    errS = geom.integrate(np.fabs(eS-eSS), dat, opts, pars, dV)

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

    for i, f in enumerate(filenames):
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

    makeErrPlot(t, nx, errRho, "magnetosonicwave_errRho", r"$L_1(\rho)$")
    makeErrPlot(t, nx, errP, "magnetosonicwave_errP", r"$L_1(P)$")
    makeErrPlot(t, nx, errVx, "magnetosonicwave_errVx", r"$L_1(v_x)$")
    makeErrPlot(t, nx, errVy, "magnetosonicwave_errVy", r"$L_1(v_y)$")
    makeErrPlot(t, nx, errVz, "magnetosonicwave_errVz", r"$L_1(v_z)$")
    makeErrPlot(t, nx, errBx, "magnetosonicwave_errBx", r"$L_1(B_x)$")
    makeErrPlot(t, nx, errBy, "magnetosonicwave_errBy", r"$L_1(B_y)$")
    makeErrPlot(t, nx, errBz, "magnetosonicwave_errBz", r"$L_1(B_z)$")
    makeErrPlot(t, nx, errS, "magnetosonicwave_errS", r"$L_1(P/\rho^\gamma)$")


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
