import sys
import h5py as h5
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discoUtil as du

def plotCheckpoint(file, vars=None):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = du.loadCheckpoint(file)
    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    piph = dat[3]

    nq = prim.shape[1]

    Zs = np.unique(z)
    z_eq = Zs[len(Zs)/2]
    eq_ind = (z==z_eq)
    title = "DISCO t = {0:.1f}".format(t)
    name = file.split('/')[-1].split('.')[0].split('_')[-1]

    if vars is None:
        vars = range(nq)

    print("   Plotting...")
    #fig, ax = plt.subplots(1,1, figsize=(10,12),subplot_kw={'projection':'polar'})
    fig, ax = plt.subplots(1,1, figsize=(4,3),subplot_kw={'projection':'polar'})
    #ax.pcolormesh(rjph, thkph, primPhi0[:,:,q], cmap=plt.cm.inferno)

    Br = prim[:,5]
    Bp = prim[:,6]
    B2 = (Br*Br + Bp*Bp) / 1.0e-8

    #vmin = B2.min()
    #vmax = B2.max()
    vmin = 0.0
    vmax = 1.05

    for i, R in enumerate(np.unique(r)[:-2]):
        ind = r==R
        phif = np.empty(len(piph[ind])+1)
        phif[1:] = piph[ind]
        phif[0] = piph[ind][-1]
        C = ax.pcolormesh(phif, rjph[i:i+2], B2[None,ind], 
                cmap=plt.cm.inferno, vmin=vmin, vmax=vmax)
    #ax.plot(np.linspace(0, 2*np.pi, 100), 2*np.ones(100), color='grey', lw=6, ls='-')
    #ax.plot(np.linspace(0, 2*np.pi, 100), 6*np.ones(100), color='grey', lw=6, ls='--')
    #ax.set_theta_direction(-1)
    #ax.set_theta_offset(0.5*np.pi)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    #ax.tick_params(axis='x', which='both', bottom='off', top='off', 
    #                labelbottom='off')
    #ax.tick_params(axis='y', which='both', bottom='off', top='off', 
    #                labelbottom='off')
    #fig.suptitle(title, fontsize=24)
    cb = fig.colorbar(C)
    cb.set_label(r"$B^2$ ($B_0^2$)")
    fig.tight_layout()
    plotname = "plot_eq_{0:s}_lin_B2.png".format(name)
    print("   Saving {0:s}...".format(plotname))
    fig.savefig(plotname, dpi=300)
    plt.close(fig)
        

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Makes plots of Disco prims as a function of r.")
        print("usage: python plotDiscoR.py <checkpoint.h5 ...>")
        sys.exit()

    files = sys.argv[1:]
    for f in files:
        plotCheckpoint(f, [5,6])
