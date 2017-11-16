from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import tracerUtil as tu
import argparse as ag


blue = (114/255.0, 158/255.0, 206/255.0)
orange = (255/255.0, 158/255.0, 74/255.0)
green = (103/255.0, 191/255.0, 92/255.0)
red = (237/255.0, 102/255.0, 93/255.0)
violet = (173/255.0, 139/255.0, 201/255.0)


def getTracerTraj( inFile, tr_ID ):

    f = open( inFile, "r" )
    pos = []

    line = "start"
    while line !="":
        line = f.readline()
        if line == "":
            break
        Ntr = int( line )
        skip = f.readline()
        for i in range( Ntr ):
            line = f.readline()
            step, time, Id, Type, xxx, yyy, zzz, rrr, phi = tu.getTracerData(line)
            if Id == tr_ID:
                pos.append( [step,time,xxx,yyy,zzz] )
        print "Checking step ", step

    f.close()
    x = np.asarray( pos  )
    filename = "tracer_{0:d}.xyz".format( tr_ID )
    np.savetxt( filename, x)
    return #x

def plotTracerTraj( tr_file, dim ):

    #filename = "tracer{0:d}.xyz".format( tr_ID )
    print("Loading {0:s}...".format(tr_file))
    x = np.loadtxt( tr_file )
    Id = tr_file.split('/')[-1].split('.')[0].split('_')[-1]

    step = x[:,0]
    time = x[:,1]
    xxx  = x[:,2]
    yyy  = x[:,3]
    zzz  = x[:,4]

    cavity = plt.Circle((0,0), 2*a, color='k', fill=False, ls='--', lw=0.5, label='Inner Cavity')
    disk   = plt.Circle((0,0), 5*a, color='k', fill=False, ls='--', lw=0.5, label='Outer Edge')

    if dim==2:
        print("   Plotting...")
        fig = plt.figure()
        ax  = fig.add_subplot(111) #, ylim=[-5,5],xlim=[-5,5])
        #ax.plot( xxx, yyy, label='tracer {0:d}'.format(tr_ID), color=blue, ls=' ', marker='o', ms=4, mew=0.5 )
        ax.plot(xxx,yyy,color=blue, ls='', marker='x')
        ax.plot( 0, 0, label='cm', color="black", marker='+', ms=8)
        ax.set_xlabel( '$x$', fontsize=18 )
        ax.set_ylabel( '$y$', fontsize=18 )
        ax.add_artist(cavity)
        ax.add_artist( disk )
        #ax.legend( loc='best', prop={'size':12} )

    elif dim==3:
        print "Still need to implement 3D tracer trajectory plotting"
    else:
        print "Dimension needs to be either 2 or 3"

    plotname = "tracer_traj_{0:s}.png".format( Id )
    print("   Saving {0:s}...".format(plotname))
    plt.savefig(plotname)
    plt.close(fig)

    return fig

if __name__ == '__main__':

    parser = ag.ArgumentParser(description="Plot Trajectories of DISCO Tracers.")
    parser.add_argument('tracers', nargs='+',
                            help="Tracer trajectory (.xyz) files to plot.")

    args = parser.parse_args()
    files = args.tracers

    a = 1
    dim = 2

    for f in files:
        plotTracerTraj( f, dim )

    #IF TRACER TRAJ FILE DOESNT EXIST: make it
    #inFile = "tracers.xyz"
    #getTracerTraj(  inFile, Id )

    #ELSE: plot
    #fig = plotTracerTraj( Id, dim )
    #plotTracerTraj( Id, dim, ax )
    #plt.show()
