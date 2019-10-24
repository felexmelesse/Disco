import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.plot as plot
import discopy.geom as geom

figsize = (9,12)
figsizeR = (8,6)

def makeGridPlots(rjph, zkph, pars, opts, title, name, fv, bounds, texlabels, 
                    labels, linear, log):

    for i in range(len(fv)):
        if linear[i]:
            fig, ax = plt.subplots(1,1, figsize=figsize)
            plot.plotPhiSlice(fig, ax, rjph, zkph, fv[i], texlabels[i], 
                            pars, opts, vmin=bounds[i,0], vmax=bounds[i,1],
                            log=False)
            fig.suptitle(title, fontsize=24)
            plotname = "plot_torus_{0:s}_lin_{1:s}.png".format(
                            name, labels[i])
            print("Saving " + plotname)
            fig.savefig(plotname, dpi=200)
            plt.close(fig)
        if log[i]:
            fig, ax = plt.subplots(1,1, figsize=figsize)
            plot.plotPhiSlice(fig, ax, rjph, zkph, fv[i], texlabels[i], 
                            pars, opts, vmin=bounds[i,0], vmax=bounds[i,1],
                            log=True)
            fig.suptitle(title, fontsize=24)
            plotname = "plot_torus_{0:s}_log_{1:s}.png".format(
                            name, labels[i])
            print("Saving " + plotname)
            fig.savefig(plotname)
            plt.close(fig)

def makeRadialPlots(r, title, name, fv, bounds, texlabels, labels, linear, 
                    log):

    for i in range(len(fv)):
        if linear[i]:
            fig, ax = plt.subplots(1,1, figsize=figsizeR)
            ax.plot(r, fv[i], '+')
            ax.set_ylim(bounds[i,0], bounds[i,1])
            ax.set_xlabel(r'$r$')
            ax.set_ylabel(texlabels[i])
            fig.suptitle(title, fontsize=24)
            plotname = "plot_torus_r_{0:s}_lin_{1:s}.png".format(
                            name, labels[i])
            print("Saving " + plotname)
            fig.savefig(plotname, dpi=200)
            plt.close(fig)
        if log[i]:
            fig, ax = plt.subplots(1,1, figsize=figsizeR)
            ax.plot(r, fv[i], '+')
            ax.set_ylim(bounds[i,0], bounds[i,1])
            ax.set_xlabel(r'$r$')
            ax.set_ylabel(texlabels[i])
            ax.set_yscale("log")
            fig.suptitle(title, fontsize=24)
            plotname = "plot_torus_r_{0:s}_log_{1:s}.png".format(
                            name, labels[i])
            print("Saving " + plotname)
            fig.savefig(plotname)
            plt.close(fig)

def makeTemporalPlots(t, fv, texlabels, labels, linear, log):

    for i in range(len(fv)):
        if linear[i]:
            fig, ax = plt.subplots(1,1, figsize=figsizeR)
            ax.plot(t, fv[i], '+')
            ax.set_xlabel(r'$t$')
            ax.set_ylabel(texlabels[i])
            plotname = "plot_torus_t_lin_{0:s}.png".format(labels[i])
            print("Saving " + plotname)
            fig.savefig(plotname, dpi=200)
            plt.close(fig)
        if log[i]:
            fig, ax = plt.subplots(1,1, figsize=figsizeR)
            ax.plot(t, fv[i], '+')
            ax.set_xlabel(r'$t$')
            ax.set_ylabel(texlabels[i])
            ax.set_yscale("log")
            plotname = "plot_torus_t_log_{0:s}.png".format(labels[i])
            print("Saving " + plotname)
            fig.savefig(plotname, dpi=200)
            plt.close(fig)


def analyze(file, bounds=None, boundsR=None, checkBounds=False, plot=True):

    print("Loading " + file)

    pars = util.loadPars(file)
    opts = util.loadOpts(file)
    t, r, phi, th, prim, dat = util.loadCheckpoint(file)
    rjph = dat[0]
    zkph = dat[1]
    iPhi0 = dat[6]

    title = "DISCO t = {0:.1f}".format(t)
    name = file.split('/')[-1].split('.')[0].split('_')[-1]

    rho = prim[:,0]
    P = prim[:,1]
    vr = prim[:,2]
    vp = prim[:,3]
    vt = prim[:,4]
    Br = prim[:,5]
    Bp = prim[:,6]
    Bt = prim[:,7]
    gam = pars['Adiabatic_Index']

    sinth = np.sin(th)

    cs = np.sqrt(gam*P/rho)
    v2 = vr*vr + r*r*vt*vt + r*r*sinth*sinth*vp*vp
    B2 = Br*Br + Bt*Bt + Bp*Bp
    Tg = P/rho
    ibeta = 0.5*B2 / P
    cA2 = B2/rho
    cA = np.sqrt(cA2)
    l = r*sinth*vp
    machS = np.sqrt(v2)/cs
    machA = np.sqrt(v2)/np.maximum(cA, 1.0e-5)
    mach = np.sqrt(v2)/np.sqrt(cs*cs + cA2)

    R = geom.getCentroid(rjph[:-1], rjph[1:], 1, opts)

    j = rho*r*r*sinth*sinth*vp
    e = 0.5*rho*v2 + P/(gam-1) + 0.5*B2

    fv = [cs[iPhi0], B2[iPhi0], Tg[iPhi0], ibeta[iPhi0], cA[iPhi0], l[iPhi0],
            mach[iPhi0], machS[iPhi0], machA[iPhi0]]
    texlabels = [r"$c_s$", r"$B^2$", r"$P/\rho$", r"$\beta^{-1}$", r"$c_A$",
                    r"$\ell$", r"$\mathcal{M}$", r"$\mathcal{M}_s$", 
                    r"$\mathcal{M}_A$"]
    labels = ["cs", "B2", "PoRho", "ibeta", "cA", "l", "mach", "machS",
                "machA"]
    linear = [True, True, True, True, True, True, True, True, True]
    log = [True, True, True, True, True, False, True, True, True]

    csmin = cs.min()
    csmax = cs.max()
    B2min = B2.min()
    if B2min < 1.0e-10:
        B2min = 1.0e-10
    B2max = B2.max()
    Tgmin = Tg.min()
    Tgmax = Tg.max()
    ibmin = ibeta.min()
    ibmax = ibeta.max()
    if ibmin < 1.0e-10:
        ibmin = 1.0e-10
    cAmin = cA.min()
    cAmax = cA.max()
    if cAmin < 1.0e-5:
        cAmin = 1.0e-5
    lmin = l.min()
    lmax = l.max()
    machmin = mach.min()
    machmax = mach.max()
    machSmin = machS.min()
    machSmax = machS.max()
    machAmin = machA.min()
    machAmax = machA.max()
    if machmin < 1.0e-5:
        machmin = 1.0e-5
    if machSmin < 1.0e-5:
        machSmin = 1.0e-5
    if machAmin < 1.0e-5:
        machAmin = 1.0e-5

    mybounds = np.array([ 
                    [csmin, csmax],
                    [B2min, B2max],
                    [Tgmin, Tgmax],
                    [ibmin, ibmax],
                    [cAmin, cAmax],
                    [lmin, lmax],
                    [machmin, machmax],
                    [machSmin, machSmax],
                    [machAmin, machAmax]])

    if bounds is None:
        bounds = mybounds
    elif checkBounds:
        bounds[:,0] = np.minimum(bounds[:,0], mybounds[:,0])
        bounds[:,1] = np.maximum(bounds[:,1], mybounds[:,1])

    if plot:
        makeGridPlots(rjph, zkph, pars, opts, title, name, fv, bounds, 
                        texlabels, labels, linear, log)

    dAr = geom.getDA1(r, dat[3], dat[1], dat, opts, pars)
    Mdot = geom.integrateTrans1(r, rho*vr, dat, opts, pars, dAr)
    phiB = geom.integrateTrans1(r, 0.5*np.abs(Br), dat, opts, pars, dAr)
    Jdot = geom.integrateTrans1(r, j*vr - r*sinth*Bp*Br, dat, opts, pars, dAr)
    Edot = geom.integrateTrans1(r, (e + P+0.5*B2)*vr, dat, opts, pars, dAr)
    Mdotmin = Mdot.min()
    Mdotmax = Mdot.max()
    Jdotmin = Jdot.min()
    Jdotmax = Jdot.max()
    phiBmin = phiB.min()
    phiBmax = phiB.max()
    Edotmin = Edot.min()
    Edotmax = Edot.max()

    fvR = [Mdot, Jdot, phiB, Edot]
    texlabelsR = [r"$\dot{M}$", r"$\dot{J}$", r"$\phi_B$", r"$\dot{E}$"]
    labelsR = ["Mdot", "Jdot", "phiB", "Edot"]
    linearR = [True, True, True, True]
    logR = [False, False, False, False]
    myboundsR = np.array([ 
                    [Mdotmin, Mdotmax],
                    [Jdotmin, Jdotmax],
                    [phiBmin, phiBmax],
                    [Edotmin, Edotmax]])

    if boundsR is None:
        boundsR = myboundsR
    elif checkBounds:
        boundsR[:,0] = np.minimum(boundsR[:,0], myboundsR[:,0])
        boundsR[:,1] = np.maximum(boundsR[:,1], myboundsR[:,1])

    if plot:
        makeRadialPlots(R, title, name, fvR, boundsR, 
                        texlabelsR, labelsR, linearR, logR)

    dV = geom.getDV(dat, opts, pars)
    M = geom.integrate(rho, dat, opts, pars, dV)
    J = geom.integrate(j, dat, opts, pars, dV)
    E = geom.integrate(e, dat, opts, pars, dV)
    EB = geom.integrate(0.5*B2, dat, opts, pars, dV)

    return t, (M, J, E, EB)
 

def analyzeAll(filenames):

    T = []
    dat = []

    bounds = np.empty((9,2))
    bounds[:,0] = np.inf
    bounds[:,1] = -np.inf
    boundsR = np.empty((4,2))
    boundsR[:,0] = np.inf
    boundsR[:,1] = -np.inf

    for f in filenames:
        analyze(f, bounds=bounds, boundsR=boundsR, checkBounds=True, 
                plot=False)

    for f in filenames:
        t, ret = analyze(f, bounds=bounds, boundsR=boundsR, 
                            checkBounds=False, plot=True)
        T.append(t)
        dat.append(ret)

    T = np.array(T)
    dat = np.array(dat).T

    texlabels = [r"$M$", r"$J$", r"$E$", r"$E_B$"]
    labels = ["M", "J", "E", "enB"]
    linear = [True, True, True, True]
    log = [True, True, True, True]

    makeTemporalPlots(T, dat, texlabels, labels, linear, log)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need some checkpoints please!")
        sys.exit()

    filenames = sys.argv[1:]

    analyzeAll(filenames)

