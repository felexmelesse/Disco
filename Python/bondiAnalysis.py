import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import discoUtil as du
import calc as ca

def analyzeSingle(filename):

    opts = du.loadOpts(filename)
    pars = du.loadPars(filename)
    print("Loading " + filename)
    t, x1, x2, x3, prim, dat = du.loadCheckpoint(filename)

    #assume spherical geom for now
    r = x1
    phi = x2
    th = x3

    GM = 1.0
    Mdot = 1.0
    rs = pars['Init_Par1']
    gam = pars['Adiabatic_Index']

    rho = prim[:,0]
    P = prim[:,1]
    vr = prim[:,2]
    vp = prim[:,3]
    vt = prim[:,4]

    cs = np.sqrt(gam*P/rho)

    us2 = GM/(2*rs)
    as2 = us2
    us = -math.sqrt(us2)
    rhos = -Mdot / (4*math.pi*rs*rs*us)
    K = as2/(gam*math.pow(rhos, gam-1))
    a02 = 0.5*(5.0-3*gam)*as2
    rho0 = math.pow(a02 / (gam*K), 1.0/(gam-1))

    R = np.logspace(math.log10(r.min()), math.log10(r.max()), num=200, base=10)

    rhoS, uS, PS = ca.bondi_newt(Mdot, GM, gam, rho0, R)
    csS = np.sqrt(gam*PS/rhoS)

    figname = "bondi_{0:s}.png".format('.'.join(
        filename.split('/')[-1].split('.')[:-1]))

    if opts['NUM_C'] == 5:

        fig, ax = plt.subplots(2,3,figsize=(14,9))
        ax[0,0].plot(r, rho, 'k+')
        ax[0,0].plot(R, rhoS)
        ax[0,0].set_xscale('log')
        ax[0,0].set_yscale('log')
        ax[0,0].set_ylabel(r'$\rho$')
        ax[0,1].plot(r, P, 'k+')
        ax[0,1].plot(R, PS)
        ax[0,1].set_xscale('log')
        ax[0,1].set_yscale('log')
        ax[0,1].set_ylabel(r'$P$')
        ax[0,2].plot(r, cs+vr, 'k+')
        ax[0,2].plot(R, csS+uS)
        ax[0,2].axhline(0.0, lw=2, ls='--', color='grey')
        ax[0,2].axvline(rs, lw=2, ls='--', color='grey')
        ax[0,2].set_xscale('log')
        ax[0,2].set_yscale('linear')
        ax[0,2].set_ylabel(r'$u^r + c_s$')
        ax[1,0].plot(r, -vr, 'k+')
        ax[1,0].plot(R, -uS)
        ax[1,0].set_xscale('log')
        ax[1,0].set_yscale('log')
        ax[1,0].set_ylabel(r'$u^r$')
        ax[1,0].set_xlabel(r'$r$')
        ax[1,1].plot(r, vp, 'k+')
        ax[1,1].plot(R, np.zeros(R.shape))
        ax[1,1].set_xscale('log')
        ax[1,1].set_yscale('linear')
        ax[1,1].set_ylabel(r'$u^\phi$')
        ax[1,1].set_xlabel(r'$r$')
        ax[1,2].plot(r, vt, 'k+')
        ax[1,2].plot(R, np.zeros(R.shape))
        ax[1,2].set_xscale('log')
        ax[1,2].set_yscale('linear')
        ax[1,2].set_ylabel(r'$u^\theta$')
        ax[1,2].set_xlabel(r'$r$')

    elif opts['NUM_C'] == 8:
        Br = prim[:,5]
        Bp = prim[:,6]
        Bt = prim[:,7]
    
        B0 = pars['Init_Par3']
        BS = B0 * rs*rs/(R*R)
        cAS = BS/np.sqrt(rhoS)
        cfS = np.sqrt(csS*csS + cAS*cAS)

        fig, ax = plt.subplots(3,3,figsize=(14,12))
        ax[0,0].plot(r, rho, 'k+')
        ax[0,0].plot(R, rhoS)
        ax[0,0].set_xscale('log')
        ax[0,0].set_yscale('log')
        ax[0,0].set_ylabel(r'$\rho$')
        ax[0,1].plot(r, P, 'k+')
        ax[0,1].plot(R, PS)
        ax[0,1].set_xscale('log')
        ax[0,1].set_yscale('log')
        ax[0,1].set_ylabel(r'$P$')
        ax[0,2].plot(r, cs+vr, 'k+')
        ax[0,2].plot(R, csS+uS)
        ax[0,2].plot(R, cfS+uS)
        ax[0,2].axhline(0.0, lw=2, ls='--', color='grey')
        ax[0,2].axvline(rs, lw=2, ls='--', color='grey')
        ax[0,2].set_xscale('log')
        ax[0,2].set_yscale('linear')
        ax[0,2].set_ylabel(r'$u^r + c_s$')
        ax[1,0].plot(r, -vr, 'k+')
        ax[1,0].plot(R, -uS)
        ax[1,0].set_xscale('log')
        ax[1,0].set_yscale('log')
        ax[1,0].set_ylabel(r'$u^r$')
        ax[1,1].plot(r, vp, 'k+')
        ax[1,1].plot(R, np.zeros(R.shape))
        ax[1,1].set_xscale('log')
        ax[1,1].set_yscale('linear')
        ax[1,1].set_ylabel(r'$u^\phi$')
        ax[1,2].plot(r, vt, 'k+')
        ax[1,2].plot(R, np.zeros(R.shape))
        ax[1,2].set_xscale('log')
        ax[1,2].set_yscale('linear')
        ax[1,2].set_ylabel(r'$u^\theta$')
        ax[2,0].plot(r, Br, 'k+')
        ax[2,0].plot(R, BS)
        ax[2,0].set_xscale('log')
        ax[2,0].set_yscale('log')
        ax[2,0].set_ylabel(r'$B^r$')
        ax[2,0].set_xlabel(r'$r$')
        ax[2,1].plot(r, Bp, 'k+')
        ax[2,1].plot(R, np.zeros(R.shape))
        ax[2,1].set_xscale('log')
        ax[2,1].set_yscale('linear')
        ax[2,1].set_ylabel(r'$B^\phi$')
        ax[2,1].set_xlabel(r'$r$')
        ax[2,2].plot(r, Bt, 'k+')
        ax[2,2].plot(R, np.zeros(R.shape))
        ax[2,2].set_xscale('log')
        ax[2,2].set_yscale('linear')
        ax[2,2].set_ylabel(r'$B^\theta$')
        ax[2,2].set_xlabel(r'$r$')

    fig.tight_layout()
    print("Saving " + figname)
    fig.savefig(figname)
    plt.close(fig)

def analyze(filenames):

    for f in filenames:
        analyzeSingle(f)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need me some checkpoints bub!")
        sys.exit()

    filenames = sys.argv[1:]
    analyze(filenames)

