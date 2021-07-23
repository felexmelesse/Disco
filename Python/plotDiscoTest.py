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
    fig, ax = plt.subplots(1,1, figsize=(12,9))
    ke = 0.5*(rho*(vr*vr+ omega*omega*r*r))
    plot.plotZSlice(fig, ax, rjph, piph, r, ke, z, "KE",
                                pars, opts)
    figname = "KE_%s.png" % name
    print("saving",figname)
    fig.savefig(figname)
    plt.close(fig)

    fig, ax = plt.subplots(1,1, figsize=(12,9))
    Angular_M = rho*(omega*r*r)
    plot.plotZSlice(fig, ax, rjph, piph, r, Angular_M, z, "Angular M",
                                pars, opts)
    figname = "Angular_M_%s.png" % name
    print("saving",figname)
    fig.savefig(figname)
    plt.close(fig)

    energy = (ke/rho) - (GM/r)
    j = Angular_M/rho

    fig, ax = plt.subplots(1,1, figsize=(12,9))
    e2 = 1 + (2*j*j*energy)/(GM*GM)
    e2[e2<0] = 0
    e = np.sqrt(e2)
    plot.plotZSlice(fig, ax, rjph, piph, r, e, z, "e",
                                pars, opts)
    figname = "e_%s.png" % name
    print("saving",figname)
    fig.savefig(figname)
    plt.close(fig)

    fig, ax = plt.subplots(1,1, figsize=(12,9))
    Div_v = geom.calculateDivV(r, phi, z, vr, omega, vz, dat, opts, pars)
    plot.plotZSlice(fig, ax, rjph, piph, r, Div_v, z, "Div_v",
                                pars, opts, vmin=-0.5, vmax=0.5)
    figname = "Div_v_%s.png" % name
    print("saving",figname)
    fig.savefig(figname)
    plt.close(fig)



if __name__=="__main__":

    for f in sys.argv[1:]:
        plotCheckpoint(Path(f))
