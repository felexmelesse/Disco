import sys
import math
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import discoUtil as du

xscale = "log"
yscale = "log"
GAM = 1.66666666667
RMIN = 0.0
RMAX = np.inf
M = 1.0


def bondiExact(Mdot, Rs, M, R):

    import calc as ca
    us2 = M / (2*Rs)
    as2 = us2 / (1-3*us2)
    us = -math.sqrt(us2)
    rhos = -Mdot / (4*np.pi*Rs*Rs*us)
    a0 = math.sqrt((GAM-1) * (1 - math.sqrt(1+3*as2) * abs(1-as2/(GAM-1))))

    rho, uSC, P = ca.bondi_rel(Mdot, M, GAM, a0, R)

    u = uSC  # Kerr-Schild u^R == Schwarzschild u^R

    return rho, u, P

def plotCheckpoint(file, plotExact=False, exactMdot=0.0, exactRs=0.0):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim = du.loadCheckpoint(file)
    if RMAX > 0:
        R2 = r*r + z*z
        ind = (R2 < RMAX*RMAX) * (R2 > RMIN*RMIN)
        r = r[ind]
        phi = phi[ind]
        z = z[ind]
        prim = prim[ind,:]

    rho = prim[:,0]
    P = prim[:,1]
    lr = prim[:,2]
    lp = prim[:,3]
    lz = prim[:,4]
    try:
        Br = prim[:,5]
        Bp = prim[:,6]/r
        Bz = prim[:,7]
    except:
        Br = np.zeros(r.shape)
        Bp = np.zeros(r.shape)
        Bz = np.zeros(r.shape)

    th = np.arctan2(r, z)
    R = np.sqrt(r*r+z*z)
    sinth = r/R
    costh = z/R

    try:
        Phi0 = prim[:,9]
        Phi1 = prim[:,10]
        Phi2 = prim[:,11]
    except IndexError:
        Phi0 = np.zeros(r.shape)
        Phi1 = np.zeros(r.shape)
        Phi2 = np.zeros(r.shape)

    print(prim.shape)

    al = 1.0 / np.sqrt(1.0 + 2*M/R)
    ber = 2*M/R * sinth / (1.0 + 2*M/R)
    bep = np.zeros(r.shape)
    bez = 2*M/R * costh / (1.0 + 2*M/R)
    gamrr = 1.0 + 2*M/R * sinth*sinth
    gamrp = np.zeros(r.shape)
    gamrz = 2.0*M/R * sinth*costh
    gampp = r*r
    gampz = np.zeros(r.shape)
    gamzz = 1.0 + 2*M/R * costh*costh
    igamrr = gamzz / (1.0 + 2*M/R)
    igamrp = np.zeros(r.shape)
    igamrz = -gamrz / (1.0 + 2*M/R)
    igampp = 1.0/(r*r)
    igampz = np.zeros(r.shape)
    igamzz = gamrr / (1.0 + 2*M/R)


    B2 = gamrr*Br*Br + 2*gamrp*Br*Bp + 2*gamrz*Br*Bz \
            + gampp*Bp*Bp + gampz*Bp*Bz + gamzz*Bz*Bz
    u2 = igamrr*lr*lr + 2*igamrp*lr*lp + 2*igamrz*lr*lz \
            + igampp*lp*lp + igampz*lp*lz + igamzz*lz*lz
    w2 = 1.0 + u2
    w = np.sqrt(w2)
    uB = lr*Br + lp*Bp + lz*Bz
    b2 = (B2 + uB*uB) / w2

    u0 = w / al
    l0 = -al*w + ber*lr + bep*lp + bez*lz
    ur = igamrr*lr + igamrp*lp + igamrz*lz - ber*u0
    up = igamrp*lr + igampp*lp + igampz*lz - bep*u0
    uz = igamrz*lr + igampz*lp + igamzz*lz - bez*u0
    uR = sinth*ur + costh*uz
    BR = sinth*Br + costh*Bz

    rhoh = rho + GAM/(GAM-1.0)*P
    s = np.log(P * np.power(rho, -GAM)) / (GAM-1.0)

    cs = np.sqrt(GAM*P/rhoh)
    cA = np.sqrt(b2/(rhoh+b2))
    Ma = np.sqrt(u2) / (cs/np.sqrt(1-cs*cs))
    Mdot = -R*R*rho*uR



    print("   Plotting...")
    nq = prim.shape[1]

    col = (z-z.min()) / (z.max() - z.min())
    col = np.log((r-r.min()+0.1) / (r.max() - r.min()))
    cm = plt.cm.viridis

    fig, ax = plt.subplots(4,4,figsize=(16,12))
    du.plotAx(ax[0,0], R, rho, xscale, yscale, r"$R$", r"$\rho$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[0,1], R, P, xscale, yscale, r"$R$", r"$P$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[0,2], R, P/rho, xscale, yscale, r"$R$", r"$P/\rho$",
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[0,3], R, cs, xscale, yscale, r"$R$", r"$c_s$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[1,0], R, uR, xscale, "linear", r"$R$", r"$u^R$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[1,1], R, up, xscale, "linear", r"$R$", r"$u^\phi$",
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[1,2], R, s, xscale, yscale, r"$R$", r"$s$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[1,3], R, Ma, xscale, yscale, r"$R$", r"$\mathcal{M}$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[2,0], R, BR, xscale, "linear", r"$R$", r"$B^R$", 
                marker='+', c=col, cmap=cm)
    ax[2,0].set_ylim(BR.min(), BR.max())
    du.plotAx(ax[2,1], R, Bp, xscale, "linear", r"$R$", r"$B^\phi$",
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[2,2], R, b2/P, xscale, yscale, r"$R$", r"$b^2/P$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[2,3], R, cA, xscale, yscale, r"$R$", r"$c_A$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[3,0], R, Mdot, xscale, "linear", r"$R$", 
                r"$\dot{M} = R^2\rho u^R$", 
                marker='+', c=col, cmap=cm)

    if plotExact:
        try:
            RR = np.logspace(math.log10(2.0*M), math.log10(R.max()), 100)
            rhoE, uE, PE = bondiExact(exactMdot, exactRs, M, RR)
            csE = np.sqrt(GAM * PE / (rhoE + GAM/(GAM-1)*PE))
            sE = np.log(PE * np.power(rhoE, -GAM)) / (GAM-1.0)
            ax[0,0].plot(RR, rhoE, 'r-')
            ax[0,1].plot(RR, PE, 'r-')
            ax[0,2].plot(RR, PE/rhoE, 'r-')
            ax[0,3].plot(RR, csE, 'r-')
            ax[1,0].plot(RR, uE, 'r-')
            ax[1,2].plot(RR, sE, 'r-')
            ax[3,0].plot(RR, exactMdot/(4*np.pi)*np.ones(RR.shape), 'r-')
        except ImportError:
            pass

    title = "DISCO t = {0:.3g}".format(t)

    #fig.suptitle(title, fontsize=18)

    plt.tight_layout()

    name = file.split('/')[-1].split('.')[0].split('_')[-1]
    plotname = "plot_bondi_{0:s}.png".format(name)
    
    print("   Saving {0:s}...".format(plotname))
    fig.savefig(plotname)

    plt.close(fig)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Makes plots of Disco prims as a function of r.")
        print("usage: python plotDiscoR.py <checkpoint.h5 ...>")
        sys.exit()

    files = sys.argv[1:]
    for f in files:
        plotCheckpoint(f, plotExact=True, exactMdot=1.0, exactRs=10.0)
