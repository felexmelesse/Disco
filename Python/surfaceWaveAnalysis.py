import sys
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom


def analyzeSingle(filename):

    print("Loading " + filename)
    name = filename.split('/')[-1].split('.')[0].split('_')[-1]

    pars = util.loadPars(filename)
    opts = util.loadOpts(filename)
    t, x, y, z, prim, dat = util.loadCheckpoint(filename)
    gam = pars['Adiabatic_Index']

    xf = dat[0]
    X = 0.5*(xf[:-1] + xf[1:])

    Nz = pars['Num_Z']

    if Nz == 1:
        sig = geom.integrate2(prim[:, 0], dat, opts, pars)
        pi = geom.integrate2(prim[:, 1], dat, opts, pars)
        vx = geom.integrate2(prim[:, 2], dat, opts, pars) / sig
        vz = np.zeros(sig.shape)
        title = r"1D --- $\Gamma$={0:.2f}".format(gam)
    else:
        sig = geom.integrateTrans1(x, prim[:, 0], dat, opts, pars)
        pi = geom.integrateTrans1(x, prim[:, 1], dat, opts, pars)
        vx = geom.integrateTrans1(x, prim[:, 2], dat, opts, pars) / sig
        vz = geom.integrateTrans1(x, prim[:, 4], dat, opts, pars) / sig
        title = r"2D --- $\Gamma$={0:.2f}".format(gam)

    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    ax[0, 0].plot(X, sig)
    ax[0, 1].plot(X, pi)
    ax[1, 0].plot(X, vx)
    ax[1, 1].plot(X, vz)
    ax[0, 0].set_ylabel(r"$\Sigma$")
    ax[0, 1].set_ylabel(r"$\Pi$")
    ax[1, 0].set_ylabel(r"$\langle v_x \rangle$")
    ax[1, 1].set_ylabel(r"$\langle v_z \rangle$")
    ax[1, 0].set_xlabel(r"$x$")
    ax[1, 1].set_xlabel(r"$x$")
    fig.suptitle(title)
    figname = "surfaceWave_{0:s}.png".format(name)
    print("Saving " + figname)
    fig.savefig(figname)
    plt.close(fig)


def analyze(filenames):

    for f in filenames:
        analyzeSingle(f)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Need checkpoints")
        sys.exit()

    filenames = sys.argv[1:]

    analyze(filenames)
