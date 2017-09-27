import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import discoUtil as du

xscale = "log"
yscale = "log"
GAM = 1.66666666667
RMAX = 0.8

def plotCheckpoint(file):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim = du.loadCheckpoint(file)

    rho = prim[:,0]
    P = prim[:,1]
    ur = prim[:,2]
    up = prim[:,3]
    uz = prim[:,4]
    Br = prim[:,5]
    Bp = prim[:,6]
    Bz = prim[:,7]

    B2 = Br*Br + Bp*Bp + Bz*Bz
    u2 = ur*ur + r*r*up*up + uz*uz
    w2 = 1.0 + u2
    uB = ur*Br + r*up*Bp + uz*Bz
    b2 = (B2 + uB*uB) / w2
    rhoh = rho + GAM/(GAM-1.0)*P

    cs = np.sqrt(GAM*P/rhoh)
    cA = np.sqrt(b2/(rhoh+b2))
    Ma = np.sqrt(u2) / (cs/np.sqrt(1-cs*cs))

    print("   Plotting...")
    nq = prim.shape[1]



    fig, ax = plt.subplots(3,4,figsize=(16,9))
    du.plotAx(ax[0,0], r, rho, xscale, yscale, r"$r$", r"$\rho$", 'k+')
    du.plotAx(ax[0,1], r, P, xscale, yscale, r"$r$", r"$P$", 'k+')
    du.plotAx(ax[0,2], r, P/rho, xscale, yscale, r"$r$", r"$P/\rho$", 'k+')
    du.plotAx(ax[0,3], r, cs, xscale, yscale, r"$r$", r"$c_s$", 'k+')
    #if nq > 8:
    #    du.plotAx(ax[0,2], r, prim[:,8], xscale, "linear", r"$r$", r"$q$", 'k+')
    du.plotAx(ax[1,0], r, ur, xscale, "linear", r"$r$", r"$u_r$", 'k+')
    du.plotAx(ax[1,1], r, up, xscale, "linear", r"$r$", r"$u_\phi$",'k+')
    du.plotAx(ax[1,2], r, uz, xscale, "linear", r"$r$", r"$u_z$", 'k+')
    du.plotAx(ax[1,3], r, Ma, xscale, yscale, r"$r$", r"$\mathcal{M}$", 'k+')
    du.plotAx(ax[2,0], r, Br, xscale, "linear", r"$r$", r"$B^r$", 'k+')
    du.plotAx(ax[2,1], r, Bp, xscale, "linear", r"$r$", r"$B^\phi$",'k+')
    #du.plotAx(ax[2,2], r, Bz, xscale, "linear", r"$r$", r"$B^z$", 'k+')
    du.plotAx(ax[2,2], r, b2/P, xscale, yscale, r"$r$", r"$b^2/P$", 'k+')
    du.plotAx(ax[2,3], r, cA, xscale, yscale, r"$r$", r"$c_A$", 'k+')

    title = "DISCO t = {0:.3g}".format(t)

    #fig.suptitle(title, fontsize=18)

    plt.tight_layout()

    name = file.split('/')[-1].split('.')[0].split('_')[-1]
    plotname = "plot_r_{0:s}.png".format(name)
    
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
        plotCheckpoint(f)
