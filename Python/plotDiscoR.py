import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import discoUtil as du

def plotCheckpoint(file):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = du.loadCheckpoint(file)

    print("   Plotting...")
    nq = prim.shape[1]

    fig, ax = plt.subplots(2,3,figsize=(14,9))
    du.plotAx(ax[0,0], r, prim[:,0], "linear", "log", r"$r$", r"$\rho$", 'k+')
    du.plotAx(ax[0,1], r, prim[:,1], "linear", "log", r"$r$", r"$P$", 'k+')
    du.plotAx(ax[1,0], r, prim[:,2], "linear", "linear", r"$r$", r"$u_r$", 'k+')
    du.plotAx(ax[1,1], r, prim[:,3], "linear", "linear", r"$r$", r"$u_\phi$",'k+')
    du.plotAx(ax[1,2], r, prim[:,4], "linear", "linear", r"$r$", r"$u_z$", 'k+')
    if nq > 5:
        du.plotAx(ax[0,2], r, prim[:,5], "linear", "linear", r"$r$", r"$q$", 'k+')

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
