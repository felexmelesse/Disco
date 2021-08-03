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
import cmasher

C = 299792458.0
G = 6.674e-11
M_solar = 1.989e30
M_BH = 1e6*M_solar
tg_sec = G*M_BH/(C*C*C)
Rg_meter = C*tg_sec
Rg_AU = 6.6846e-12*Rg_meter
t_disk = 1e4

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

    r = r*Rg_AU
    rjph = rjph*Rg_AU
    name = (file.stem).split("_")[-1]
    fig, ax = plt.subplots(1,1, figsize=(8,6))

    plot.plotZSlice(fig, ax, rjph, piph+np.pi/2, r, rho, z, r"$\Sigma$",
                                pars, opts, log=True, cmap = plt.get_cmap('cmr.ember'))
    
    ax.set_xlabel("X (AU)")
    ax.set_ylabel("Y (AU)")
    figname = "Density_%s.png" % name
    print("saving",figname)
    fig.savefig(figname,dpi=400)
    plt.close(fig)

if __name__=="__main__":

    for f in sys.argv[1:]:
        plotCheckpoint(Path(f))
