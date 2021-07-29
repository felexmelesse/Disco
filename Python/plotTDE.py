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

def analyzecheCkpoint(file):
    
    print("Loading {0:s}...".format(str(file)))

    t, r, phi, z, prim, dat = util.loadCheckpoint(file)
    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    piph = dat[3]
    pars = util.loadPars(file)
    opts = util.loadOpts(file)

    rho = prim[:,0]
    P = prim[:,1]
    vr = prim[:,2]
    omega = prim[:,3]
    vz = prim[:,4]
    GM = 1
    R = np.unique(r)

    name = (file.stem).split("_")[-1]
 
    ke = 0.5*(rho*(vr*vr+ omega*omega*r*r))
    
    J = rho*(omega*r*r)

    sq6 = np.sqrt(6)
    Alpha = -4*(2.0+sq6)/3
    Rx = GM*(4.0*sq6 - 9)
    Ry = -4*GM*(2*sq6 - 3.0)/3.0
    Phi = Alpha*GM/r + (1-Alpha)*GM/(r-Rx) + GM*Ry/(r*r)
    energy = (ke/rho)- Phi

#    energy = (ke/rho) - (GM/r)

    j = J/rho
   
    e2 = 1 + (2*j*j*energy)/(GM*GM)
    e2[e2<0] = 0
    e = np.sqrt(e2)

    HovrR = np.sqrt(P*r/(rho*GM))


# ecc_r = np.empty(Nr)
    
    rhosum = geom.integrate2(rho, dat, opts, pars)
    Jsum = geom.integrate2(J, dat, opts, pars)
    energysum = geom.integrate2(rho*energy, dat, opts, pars)
    rho_e_sum = geom.integrate2(rho*e, dat, opts, pars)
    HovrR_sum = geom.integrate2(rho*HovrR, dat, opts, pars)

    ecc_r = rho_e_sum/rhosum
    rho_r = rhosum/(2*np.pi)
    J_r = Jsum/rhosum
    energy_r = energysum/rhosum
    HovrR_r = HovrR_sum/rhosum


    return t, R, ecc_r, rho_r, HovrR_r  

if __name__=="__main__":

    fig_e, ax_e = plt.subplots(1,1, figsize=(12,9))
    fig_rho, ax_rho = plt.subplots(1,1, figsize=(12,9))
    fig_H, ax_H = plt.subplots(1,1, figsize=(12,9))
    cmap = plt.get_cmap("plasma")

    for i, f in enumerate(sys.argv[1:]):

        t, R, ecc_r,rho_r, HovrR_r= analyzecheCkpoint(Path(f))
        label = "time = %.0f" % t
        color = cmap(i/(len(sys.argv[1:])-1))
        ax_e.plot(R,ecc_r,label = label, color = color)
        ax_rho.plot(R,rho_r,label = label, color = color)
        ax_H.plot(R,HovrR_r,label = label, color = color)

    ax_e.legend()
    ax_e.set_ylim(0,1.5)
    ax_e.set_xlabel("R")
    ax_e.set_ylabel("e")
    figname = "r_vs_e_time.png"
    print("saving",figname)
    fig_e.savefig(figname)
    plt.close(fig_e)

    ax_rho.legend()
    ax_rho.set_xlabel("R")
    ax_rho.set_ylabel(r"$\rho$")
    figname = "r_vs_rho_time.png"
    print("saving",figname)
    fig_rho.savefig(figname)
    plt.close(fig_rho)

    ax_H.legend()
    ax_H.set_ylim(0,2)
    ax_H.set_xlabel("R")
    ax_H.set_ylabel("H/R")
    figname = "r_vs_HovrR_time.png"
    print("saving",figname)
    fig_H.savefig(figname)
    plt.close(fig_H)















