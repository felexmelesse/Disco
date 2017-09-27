import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import discoUtil as du
import math
import scipy.integrate as int

blue = (31.0/255, 119.0/255, 180.0/255)
orange = (255.0/255, 127.0/255, 14.0/255)

def calc_infall(r, t, rho0, P0, gam, M):

    def u0(x):
        return 2*(math.sqrt(2.)*x - math.pow(6./x-1,1.5)) / (3*(x-2.))

    def ur(x):
        return -math.pow(6./x-1, 1.5) / 3.0
    
    def v(x):
        if x < 6.:
            return ur(x) / u0(x)
        return 0.0

    def f(x):
        return 1./v(x)

    def F(x0, x):
        return int.quad(f, x0, x)[0]

    def calc_x0(x, t, TOL=1.0e-10):
        #print("x, t = ({0:g}, {1:g})".format(x, t))
        if x >= 6.:
            return x
        a = 0.0
        b = 6.0
        x0 = 0.5*(x+b)
        dx0 = np.inf
        while abs(dx0) > TOL:
            #print("   x0 = {0:g}".format(x0))
            dx0 = -(F(x0,x)-t) / (-f(x0))
            if x0 + dx0 >= b:
                x0 = 0.5*(x0+b)
            elif x0 + dx0 <= a:
                x0 = 0.5*(x0+a)
            else:
                x0 = x0 + dx0
        return x0

    rho = np.empty(r.shape)
    P = np.empty(r.shape)
    for i,r in enumerate(r):
        if r < 6*M:
            x0 = calc_x0(r/M, t/M)
            rho[i] = rho0 * x0 * u0(x0) * v(x0) / (r/M * u0(r/M) * v(r/M))
        else:
            rho[i] = rho0

    P = P0 * np.power(rho/rho0, gam)

    return rho, P


def plotCheckpoint(file):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = du.loadCheckpoint(file)

    M = 1.0
    gam = 1.66666666667

    rho = prim[8:,0]
    P = prim[8:,1]
    lr = prim[8:,2]
    lp = prim[8:,3]
    R = r[8:]

    cs2 = gam * P / (rho + gam*P/(gam-1.))
    ucs2 = cs2 / (1-cs2)

    u2 = lr*lr / (1.0+2*M/R) + lp*lp / (R*R)
    v2 = u2 / (1.0 + u2)

    mach = np.sqrt(u2 / ucs2)
    #mach = np.sqrt(v2 / cs2)

    print("   Plotting...")
    nq = prim.shape[1]

    rho_e, P_e = calc_infall(R, t, rho[-1], P[-1], gam, M)

    fig, ax = plt.subplots(1,1,figsize=(8,6))

    l1, = ax.plot(R, rho, lw=4, color=blue, label=r'$\rho$')
    ax.plot(R, rho_e, lw=6, color='grey', alpha=0.5)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_ylim(1.0e-5, 1.0e1)
    
    ax.set_xlabel(r'$r$ $(M)$', fontsize=18)
    ax.set_ylabel(r'$\rho$', fontsize=18)
    ax.tick_params(labelsize=18)
    #ax.legend(fontsize=18)

    ax2 = ax.twinx()
    l2, = ax2.plot(R, P, lw=6, color=orange, label=r'$P$', alpha=0.5)
    ax2.plot(R, P_e, lw=6, color='grey', ls='--', alpha=0.5)
    ax2.set_yscale('log')
    ax2.set_ylim(1.0e-9, 1.0e-2)
    ax2.set_ylabel(r'$P$', fontsize=18)
    ax2.tick_params(labelsize=18)
    ax2.axvspan(1.0, 2.0, color='grey', alpha=0.5)
    fig.legend([l1, l2], [r'$\rho$', r'$P$'], fontsize=18, loc=(0.7,0.2))

    #ax2 = ax.twinx()
    #l2, = ax2.plot(R, mach, lw=4, color=orange, label=r'$\mathcal{M}$')
    #ax2.set_yscale('log')
    #ax2.set_ylim(1.0, 1.0e3)
    #ax2.set_ylabel(r'$\mathcal{M}$', fontsize=18)
    #ax2.tick_params(labelsize=18)
    #fig.legend([l1, l2], [r'$\rho$', r'$\mathcal{M}$'], fontsize=18)

    #du.plotAx(ax[0,0], r, prim[:,0], "log", "log", r"$r$", r"$\rho$", 'k+')
    #du.plotAx(ax[0,1], r, prim[:,1], "log", "log", r"$r$", r"$P$", 'k+')
    #du.plotAx(ax[1,0], r, prim[:,2], "linear", "linear", r"$r$", r"$u_r$", 'k+')
    #du.plotAx(ax[1,1], r, prim[:,3], "linear", "linear", r"$r$", r"$u_\phi$",'k+')
    #du.plotAx(ax[1,2], r, prim[:,4], "linear", "linear", r"$r$", r"$u_z$", 'k+')
    #if nq > 5:
    #    du.plotAx(ax[0,2], r, prim[:,5], "linear", "linear", r"$r$", r"$q$", 'k+')

    title = "DISCO t = {0:.1f}".format(t)

    fig.suptitle(title, fontsize=18)

    plt.tight_layout()

    name = file.split('/')[-1].split('.')[0].split('_')[-1]
    plotname = "plot_rhomach_{0:s}.png".format(name)
    
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
