import sys
import h5py as h5
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discoUtil as du

def plotCheckpoint(file):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, theta, prim, dat = du.loadCheckpoint(file)
    rjph = dat[0]
    thkph = dat[1]
    primPhi0 = dat[2]

    nq = prim.shape[1]

    Thetas = np.unique(theta)
    Theta_eq = Thetas[len(Thetas)/2]
    eq_ind = (theta==Theta_eq)
    title = "DISCO t = {0:.3g}".format(t)
    name = file.split('/')[-1].split('.')[0].split('_')[-1]

    for q in range(nq):
        print("   Plotting...")
        fig, ax = plt.subplots(1,1, subplot_kw={'projection':'polar'})
        #ax.pcolormesh(rjph, thkph, primPhi0[:,:,q], cmap=plt.cm.inferno)
        C = ax.pcolormesh(thkph, rjph, primPhi0[:,:,q].T, cmap=plt.cm.inferno)
        ax.set_theta_direction(-1)
        ax.set_theta_offset(0.5*np.pi)
        ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                        labelbottom='off')
        fig.colorbar(C)
        fig.suptitle(title, fontsize=18)
        fig.tight_layout()
        plotname = "plot_Phi0_{0:s}_lin_{1:d}.png".format(name, q)
        print("   Saving {0:s}...".format(plotname))
        fig.savefig(plotname)
        plt.close(fig)
        
        if q in [0,1]:
            fig, ax = plt.subplots(1,1, subplot_kw={'projection':'polar'})
            #ax.pcolormesh(rjph, thkph, primPhi0[:,:,q], cmap=plt.cm.inferno)
            C = ax.pcolormesh(thkph, rjph, primPhi0[:,:,q].T, 
                                cmap=plt.cm.inferno, norm=mpl.colors.LogNorm())
            ax.set_theta_direction(-1)
            ax.set_theta_offset(0.5*np.pi)
            ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                            labelbottom='off')
            fig.colorbar(C)
            fig.suptitle(title, fontsize=18)
            fig.tight_layout()
            plotname = "plot_Phi0_{0:s}_log_{1:d}.png".format(name, q)
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
