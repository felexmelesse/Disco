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

    t, r, phi, th, prim = du.loadCheckpoint(file)
    if RMAX > 0:
        R2 = r*r
        ind = (R2 < RMAX*RMAX) * (R2 > RMIN*RMIN)
        r = r[ind]
        phi = phi[ind]
        th = th[ind]
        prim = prim[ind,:]

    rho = prim[:,0]
    P = prim[:,1]
    lr = prim[:,2]
    lp = prim[:,3]
    lt = prim[:,4]
    try:
        Br = prim[:,5]
        Bp = prim[:,6]/r
        Bt = prim[:,7]
    except:
        Br = np.zeros(r.shape)
        Bp = np.zeros(r.shape)
        Bt = np.zeros(r.shape)


    try:
        Phi0 = prim[:,9]
        Phi1 = prim[:,10]
        Phi2 = prim[:,11]
    except IndexError:
        Phi0 = np.zeros(r.shape)
        Phi1 = np.zeros(r.shape)
        Phi2 = np.zeros(r.shape)

    print(prim.shape)

    st = np.sin(th)
    ct = np.cos(th)

    al = 1.0 / np.sqrt(1.0 + 2*M/r)
    ber = 2*M/r / (1.0 + 2*M/r)
    bep = np.zeros(r.shape)
    beth = np.zeros(r.shape)
    gamrr = 1.0 + 2*M/r
    gamrp = np.zeros(r.shape)
    gamrth = np.zeros(r.shape)
    gampp = r*r*st*st
    gampth = np.zeros(r.shape)
    gamthth = r*r
    igamrr = 1.0 / (1.0 + 2*M/r)
    igamrp = np.zeros(r.shape)
    igamrth = np.zeros(r.shape)
    igampp = 1.0/(r*r*st*st)
    igampth = np.zeros(r.shape)
    igamthth = 1.0 / (r*r)


    B2 = gamrr*Br*Br + 2*gamrp*Br*Bp + 2*gamrth*Br*Bt \
            + gampp*Bp*Bp + gampth*Bp*Bt + gamthth*Bt*Bt
    u2 = igamrr*lr*lr + 2*igamrp*lr*lp + 2*igamrth*lr*lt \
            + igampp*lp*lp + igampth*lp*lt + igamthth*lt*lt
    w2 = 1.0 + u2
    w = np.sqrt(w2)
    uB = lr*Br + lp*Bp + lt*Bt
    b2 = (B2 + uB*uB) / w2

    u0 = w / al
    l0 = -al*w + ber*lr + bep*lp + beth*lt
    ur = igamrr*lr + igamrp*lp + igamrth*lt - ber*u0
    up = igamrp*lr + igampp*lp + igampth*lt - bep*u0
    uz = igamrth*lr + igampth*lp + igamthth*lt - beth*u0

    rhoh = rho + GAM/(GAM-1.0)*P
    s = np.log(P * np.power(rho, -GAM)) / (GAM-1.0)

    cs = np.sqrt(GAM*P/rhoh)
    cA = np.sqrt(b2/(rhoh+b2))
    Ma = np.sqrt(u2) / (cs/np.sqrt(1-cs*cs))
    Mdot = -r*r*rho*ur



    print("   Plotting...")
    nq = prim.shape[1]

    col = (th-th.min()) / (th.max() - th.min())
    #col = np.log((r-r.min()+0.1) / (r.max() - r.min()))
    cm = plt.cm.viridis

    fig, ax = plt.subplots(4,4,figsize=(16,12))
    du.plotAx(ax[0,0], r, rho, xscale, yscale, r"$R$", r"$\rho$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[0,1], r, P, xscale, yscale, r"$R$", r"$P$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[0,2], r, P/rho, xscale, yscale, r"$R$", r"$P/\rho$",
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[0,3], r, cs, xscale, yscale, r"$R$", r"$c_s$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[1,0], r, ur, xscale, "linear", r"$R$", r"$u^R$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[1,1], r, up, xscale, "linear", r"$R$", r"$u^\phi$",
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[1,2], r, s, xscale, yscale, r"$R$", r"$s$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[1,3], r, Ma, xscale, yscale, r"$R$", r"$\mathcal{M}$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[2,0], r, Br, xscale, "linear", r"$R$", r"$B^R$", 
                marker='+', c=col, cmap=cm)
    ax[2,0].set_ylim(Br.min(), Br.max())
    du.plotAx(ax[2,1], r, Bp, xscale, "linear", r"$R$", r"$B^\phi$",
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[2,2], r, b2/P, xscale, yscale, r"$R$", r"$b^2/P$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[2,3], r, cA, xscale, yscale, r"$R$", r"$c_A$", 
                marker='+', c=col, cmap=cm)
    du.plotAx(ax[3,0], r, Mdot, xscale, "linear", r"$R$", 
                r"$\dot{M} = R^2\rho u^R$", 
                marker='+', c=col, cmap=cm)

    if plotExact:
        try:
            RR = np.logspace(math.log10(2.0*M), math.log10(r.max()), 100)
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
