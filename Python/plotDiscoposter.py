import os
import sys
from pathlib import Path 
import argparse as ag
import h5py as h5
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.plot as plot
import discopy.geom as geom

def plotCheckpoint(file, vars=None, logvars=None, noGhost=False, om=None,
                    bounds=None, rmax=None, planets=False, k=None):
    
    print("Loading {0:s}...".format(str(file)))

    t, r, phi, z, prim, dat = util.loadCheckpoint(file)
    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    piph = dat[3]
    pars = util.loadPars(file)
    opts = util.loadOpts(file)

    rho = prim[:,0]
    vr = prim[:,2]
    omega = prim[:,3]
    vz = prim[:,4]
    GM = 1

    name = (file.stem).split("_")[-1]
    fig, ax = plt.subplots(1,1, figsize=(8,6))

    plot.plotZSlice(fig, ax, rjph, piph+np.pi, r, rho, z, r"$\rho$",
                                pars, opts,log=True)
    figname = "Density_%s.png" % name
    print("saving",figname)
    fig.savefig(figname)
    plt.close(fig)

if __name__=="__main__":

    for f in sys.argv[1:]:
        plotCheckpoint(Path(f))
