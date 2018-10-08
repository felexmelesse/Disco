import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom

def analyzeSingle(filename):

    opts = util.loadOpts(filename)
    pars = util.loadPars(filename)
    print("Loading " + filename)
    t, x1, x2, x3, prim, dat = util.loadCheckpoint(filename)

    B1 = prim[:,5]
    B2 = prim[:,6]
    B3 = prim[:,7]

    b2 = B1*B1 + B2*B2 + B3*B3

    eB = 0.5 * geom.integrate(b2, dat, opts, pars)

    return t, eB

def analyze(filenames):

    N = len(filenames)

    eB = np.empty(N)
    t = np.empty(N)

    for i,f in enumerate(filenames):
        t[i], eB[i] = analyzeSingle(f)

    figname = "floopMagneticEnergy.png"
    fig, ax = plt.subplots(1,1)
    ax.plot(t, eB/eB[0], 'k+')
    ax.set_xlabel("t")
    ax.set_ylabel(r"$e_B(t)/e_B(t=0)$")
    ax.set_ylim(0.85,1.0)
    print("Saving " + figname)
    fig.savefig("floopMagneticEnergy.png")

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need a checkpoint dude!")
        sys.exit()

    filenames = sys.argv[1:]
    analyze(filenames)
