import sys
import numpy as np
import discopy.util as util
import discopy.geom as geom

def calc(filename, makePlot=False, trim=False):

    t, r, phi, z,prim, dat = util.loadCheckpoint(filename)
    opts = util.loadOpts(filename)
    pars = util.loadPars(filename)

    divB = geom.calcDivB(dat, opts, pars)

    IdivB = geom.integrate(np.abs(divB), dat, opts, pars)

    print(IdivB)

    if makePlot:
        import discopy.plot as plot
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1,1, figsize=(9,9))

        rjph = dat[0]
        zkph = dat[1]
        idPhi0 = dat[6]
        index = dat[8]
        nr = len(rjph)-1
        nz = len(zkph)-1
        divB0 = np.empty((nz, nr))
        for k in range(nz):
            for j in range(nr):
                divB0[k,j] = divB[index[k,j]]

        if trim:
            rjph = rjph[2:-2]
            zkph = zkph[2:-2]
            divB0 = divB0[2:-2,2:-2]

        plot.plotPhiSlice(fig, ax, rjph, zkph, divB0, 
                            r"$\nabla \cdot B$",
                            pars, opts)
        name = ".".join(filename.split("/")[-1].split(".")[:-1])
        figname = "plot_divB_rz_" + name + ".png"
        print("Saving " + figname)
        fig.savefig(figname)
        plt.close(fig)
        return IdivB, divB0

    return IdivB, None

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need some checkpoints!")
        sys.exit()

    makePlot = False
    if "plot" in sys.argv:
        print("Gonna plot")
        sys.argv.remove("plot")
        makePlot = True
    trim = False
    if "trim" in sys.argv:
        print("Gonna trim")
        sys.argv.remove("trim")
        trim = True
    diff = False
    if "diff" in sys.argv:
        print("Gonna diff")
        sys.argv.remove("diff")
        diff = True

    val0, divB0 = calc(sys.argv[1], makePlot, trim)
    
    for f in sys.argv[2:-1]:
        calc(f, makePlot, trim)

    val1, divB1 = calc(sys.argv[-1], makePlot, trim)

    if val0 != 0.0:
        err = (val1-val0)/val0
    else:
        err = val1

    print(err)

    if makePlot:
        diff = divB1-divB0
        import discopy.plot as plot
        import matplotlib.pyplot as plt
        
        t, r, phi, z,prim, dat = util.loadCheckpoint(sys.argv[1])
        opts = util.loadOpts(sys.argv[1])
        pars = util.loadPars(sys.argv[1])

        fig, ax = plt.subplots(1,1, figsize=(9,9))

        rjph = dat[0]
        zkph = dat[1]
        idPhi0 = dat[6]
        index = dat[8]
        nr = len(rjph)-1
        nz = len(zkph)-1

        if trim:
            rjph = rjph[2:-2]
            zkph = zkph[2:-2]

        plot.plotPhiSlice(fig, ax, rjph, zkph, diff, 
                            r"$\delta \nabla \cdot B$",
                            pars, opts)
        name = ".".join(sys.argv[1].split("/")[-1].split(".")[:-1])
        figname = "plot_divB_diff_rz_" + name + ".png"
        print("Saving " + figname)
        fig.savefig(figname)
        plt.close(fig)
