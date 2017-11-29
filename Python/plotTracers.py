from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py as h5
import tracerUtil as tu
import discoUtil as du
import argparse as ag


def plotTracers( file ):

    print("Loading {0:s}...".format(file))

    dat = du.loadCheckpoint(file)[-1]
    rpz = tu.loadTrH5Data( file )[0]
    r   = rpz[:,0]
    phi = rpz[:,1]
    z   = rpz[:,2]

    x = r*np.cos(phi)
    y = r*np.sin(phi)
    t = 0.0


    rmax = ( h5.File(file, 'r') )['Pars']['R_Max'][0]
    lim = rmax + 0.5
    disk = plt.Circle( (0,0), rmax, color='k', fill=False, ls='-', lw=0.5, label='Outer Edge')

    title = "t = {0:.1f}".format(t)
    name = file.split('/')[-1].split('.')[0].split('_')[-1]

    print("   Plotting...")
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.plot( x, y ,color='b', ls='', marker='o', mew=1, mec='w', ms=5)
    ax.set_xlabel( '$x$', fontsize=18 )
    ax.set_ylabel( '$y$', fontsize=18 )
    ax.add_artist( disk )

    planets = dat[4]
    rpl = planets[:,3]
    ppl = planets[:,4]
    xpl = rpl * np.cos(ppl)
    ypl = rpl * np.sin(ppl)
    ax.plot(xpl, ypl, color='black', ls='', marker='o', mew=0, ms=8)

    ax.set_aspect('equal')
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    plotname = "tracers_{0:s}.png".format(name)
    print("   Saving {0:s}...".format(plotname))
    fig.savefig(plotname)
    plt.close(fig)

if __name__ == '__main__':

    parser = ag.ArgumentParser(description="Plot DISCO Tracers.")
    parser.add_argument('checkpoints', nargs='+',
                            help="Checkpoint (.h5) files to plot.")

    args = parser.parse_args()
    files = args.checkpoints

    for f in files:
        plotTracers( f )

    #plt.show()
