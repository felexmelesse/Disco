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
    R = np.unique(r)

    name = (file.stem).split("_")[-1]
 
    ke = 0.5*(rho*(vr*vr+ omega*omega*r*r))
    
    J = rho*(omega*r*r)

    energy = (ke/rho) - (GM/r)
    j = J/rho
   
    e2 = 1 + (2*j*j*energy)/(GM*GM)
    e2[e2<0] = 0
    e = np.sqrt(e2)

# ecc_r = np.empty(Nr)
    
    rhosum = geom.integrate2(rho, dat, opts, pars)
    Jsum = geom.integrate2(J, dat, opts, pars)
    energysum = geom.integrate2(rho*energy, dat, opts, pars)
    rho_e_sum = geom.integrate2(rho*e, dat, opts, pars)

    ecc_r = rho_e_sum/rhosum
    rho_r = rhosum/(2*np.pi)
    J_r = Jsum/rhosum
    energy_r = energysum/rhosum

    fig, ax = plt.subplots(1,1, figsize=(12,9))
    ax.plot(R,rho_r)
    ax.set_xlabel("R")
    ax.set_ylabel(r"$\rho$")
    figname = "r_vs_rho%s.png" % name
    print("saving",figname)
    fig.savefig(figname)
    plt.close(fig)

#    for j, rr in enumerate(R):
#        ind = (rr == r)
#        sigsum = np.sum(dphi[ind]*rho[ind])
#        ecc_r[j] = np.sum(e[ind]*rho[ind]*dphi[ind])/(sigsum)

    fig, ax = plt.subplots(1,1, figsize=(12,9))
    ax.plot(R,ecc_r)
    ax.set_xlabel("R")
    ax.set_ylabel("e")
    ax.set_ylim(0,1)
    figname = "r_vs_e%s.png" % name
    print("saving",figname)
    fig.savefig(figname)
    plt.close(fig)

    fig, ax = plt.subplots(1,1, figsize=(12,9))
    ax.plot(R,J_r)
    ax.set_xlabel("R")
    ax.set_ylabel("J")
    figname = "r_vs_J%s.png" % name
    print("saving",figname)
    fig.savefig(figname)
    plt.close(fig)

    fig, ax = plt.subplots(1,1, figsize=(12,9))
    ax.plot(R,energy_r)
    ax.set_xlabel("R")
    ax.set_ylabel("Energy")
    figname = "r_vs_Energy%s.png" % name
    print("saving",figname)
    fig.savefig(figname)
    plt.close(fig)


    
#   Div_v = geom.calculateDivV(r, phi, z, vr, omega, vz, dat, opts, pars)
    
    Nr = len(R)
    jplot = [3*Nr//4, 7*Nr//8, 15*Nr//16,31*Nr//32]
   
    fig_rho, ax_rho = plt.subplots(1,1, figsize=(12,9))
    fig_e, ax_e = plt.subplots(1,1, figsize=(12,9))
    fig_J, ax_J= plt.subplots(1,1, figsize=(12,9))
    fig_energy, ax_energy = plt.subplots(1,1, figsize=(12,9))
   
    for j in jplot:
        ind = (R[j] == r)
        phi_ring = phi[ind]
        phi_ring[phi_ring > np.pi] -= 2*np.pi
        rho_ring = rho[ind]
        e_ring = e[ind]
        J_ring = (J/rho)[ind]
        energy_ring = energy[ind]
        label = "r = %0.1f" % R[j]
        ax_rho.plot(phi_ring,rho_ring,ls="",marker="+",label=label)
        ax_e.plot(phi_ring,e_ring,ls="",marker="+",label=label)
        ax_J.plot(phi_ring,J_ring,ls="",marker="+",label=label)
        ax_energy.plot(phi_ring,energy_ring,ls="",marker="+",label=label)

    ax_rho.set_xlabel(r"$\phi$")
    ax_rho.set_ylabel(r"$\rho$")

    ax_e.set_xlabel(r"$\phi$")
    ax_e.set_ylabel("e")

    ax_J.set_xlabel(r"$\phi$")
    ax_J.set_ylabel("J")

    ax_energy.set_xlabel(r"$\phi$")
    ax_energy.set_ylabel("energy")
    
    ax_rho.legend()   
    ax_e.legend() 
    ax_J.legend()
    ax_energy.legend()

    figname = "phi_vs_rho%s.png" % name
    print("saving",figname)
    fig_rho.savefig(figname)
    plt.close(fig_rho)

    figname = "phi_vs_e%s.png" % name
    print("saving",figname)
    fig_e.savefig(figname)
    plt.close(fig_e)

    figname = "phi_vs_J%s.png" % name
    print("saving",figname)
    fig_J.savefig(figname)
    plt.close(fig_J)

    figname = "phi_vs_energy%s.png" % name
    print("saving",figname)
    fig_energy.savefig(figname)
    plt.close(fig_energy)





if __name__=="__main__":

    for f in sys.argv[1:]:
        plotCheckpoint(Path(f))
