import sys
import math
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom
import discopy.plot as plot


def analyzeSingle(filename):

    opts = util.loadOpts(filename)
    pars = util.loadPars(filename)
    print("Loading " + filename)
    t, x1, x2, x3, prim, dat = util.loadCheckpoint(filename)
    name = filename.split('/')[-1].split('.')[0].split('_')[-1]

    R = np.unique(x1)

    rjph = dat[0]
    piph = dat[3]

    nx = pars['Num_R']
    gam = pars['Adiabatic_Index']

    rho = prim[:, 0]
    P = prim[:, 1]
    v1 = prim[:, 2]
    v2 = prim[:, 3]
    v3 = prim[:, 4]
    q = prim[:, 5]

    x, y, z = geom.getXYZ(x1, x2, x3, opts, pars)
    vx, vy, vz = geom.getVXYZ(x1, x2, x3, v1, v2, v3, opts)

    print("Plotting") 
    fig, ax = plt.subplots(1, 1)

    dphi = 2.0*(piph-x2)

    #for j, rr in enumerate(R):
    #  ind = (rr == x1)
    #  avgRho = np.sum(rho[ind]*dphi[ind])/np.sum(ind)
    #  rho[ind] = (rho[ind]-avgRho)/avgRho

    #plot.plotZSlice(fig, ax, rjph, piph, x1, rho, x3.mean(), r"$\rho/\langle\rho\rangle_\phi$",
    #                pars, opts)
    plot.plotZSlice(fig, ax, rjph, piph, x1, rho, x3.mean(), r"$\rho$",
                    pars, opts, vmax = 1.2)
    ax.set_aspect('equal')
    print("Saving") 
    #fig.savefig('bl_deltaRho.png')
    plotname = "plot_BL_{0:s}_{1:s}.png".format(name, 'rho')
    fig.savefig(plotname, dpi = 400)
    plt.close(fig)

    print("Plotting vr")
    fig, ax = plt.subplots(1, 1)

    vrp = v1 * x1 * np.sqrt(rho)

    v1ext = np.max([-1*vrp.min(), vrp.max()])
    slt = 0.001 * v1ext

    #plot.plotZSlice(fig, ax, rjph, piph, x1, v1, x3.mean(), r"$v_r$",
    #                pars, opts, symlog=True, symlthresh = slt, vmax = v1ext, vmin = -1*v1ext, cmap=plt.get_cmap('PRGn') )

    plot.plotZSlice(fig, ax, rjph, piph, x1, vrp, x3.mean(), r"$r\sqrt{(\Sigma)}v_r$",
                    pars, opts, symlog=True, symlthresh=slt, vmax = v1ext, vmin = -1*v1ext, cmap=plt.get_cmap('PRGn') )

    ax.set_aspect('equal')
    print("Saving") 
    plotname = "plot_BL_{0:s}_{1:s}.png".format(name, 'vr')
    fig.savefig(plotname, dpi = 400)
    #fig.savefig('bl_rVr.png')
    plt.close(fig)


    print("Calculating div")
    divV = geom.calculateCurlV(x1, x2, x3, v1, v2, v3, dat, opts, pars)

    interior = (x1 > pars['R_Min']) & (x1 < pars['R_Max'])
    minDiv = divV[interior].min()
    maxDiv = divV[interior].max()
    if np.abs(maxDiv) < np.abs(minDiv): maxDiv = -1.0*minDiv
    else : minDiv = -1 * maxDiv

    print("Plotting") 
    fig, ax = plt.subplots(1, 1)
    plot.plotZSlice(fig, ax, rjph, piph, x1, divV, x3.mean(),
                    r"$(\nabla\times v)_z$", pars, opts, symlog=True, symlthresh=0.001*maxDiv, vmin=minDiv, vmax=maxDiv, cmap='RdBu')
    ax.set_aspect('equal')
    print("Saving") 
    plotname = "plot_BL_{0:s}_{1:s}.png".format(name, 'vorticity')
    fig.savefig(plotname, dpi = 400)
    #fig.savefig('bl_vorticity.png', dpi = 1200)
    plt.close(fig)



    return t, nx

def analyze(filenames):

    N = len(filenames)

    t = np.empty(N)
    nx = np.empty(N)

    for i, f in enumerate(filenames):
        dat = analyzeSingle(f)
        t[i] = dat[0]
        nx[i] = dat[1]


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need a checkpoint dude!")
        sys.exit()

    filenames = sys.argv[1:]
    analyze(filenames)
