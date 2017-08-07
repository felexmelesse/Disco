import os
import sys
import math
import argparse as ap
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discoUtil as du
import plotDiscoEq as pde

def plotGradients(filename):

    print("  Loading "+filename)
    t, r, phi, z, prim, dat = du.loadCheckpoint(filename)
    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    piph = dat[3]
    nphi = dat[7]
    pars = du.loadPars(filename)

    k = int(zkph.shape[0]/2-1)

    zind = (z>zkph[k]) * (z<zkph[k+1])
    r = r[zind]
    z = z[zind]
    phi = phi[zind]
    prim = prim[zind,:]
    primPhi0 = primPhi0[k,:]
    piph = piph[zind]
    nphi = nphi[k,:]
    id0 = np.zeros(nphi.shape, dtype=np.int)
    id0[1:] = np.cumsum(nphi)[:-1]

    grid = (rjph, piph, id0, nphi)

    planetDat = dat[4]

    varnames, vartex, num_c, num_n = du.getVarNames(filename)
    nq = num_c + num_n

    title = "DISCO t = {0:.1f}".format(t)
    name = filename.split('/')[-1].split('.')[0].split('_')[-1]

    print("    Calc Grad ")
    PLM = pars['PLM']
    PHI_MAX = pars['Phi_Max']
    gam = pars['Adiabatic_Index']
    gradr, gradp = calcCellGrad(grid, prim, PLM, PHI_MAX)

    rho = prim[:,0]
    P = prim[:,1]
    vr = prim[:,2]
    om = prim[:,3]
    s = np.log(P * np.power(rho, -gam)) / (gam-1.0)
    s -= s.min()

    divv = gradr[:,2] + vr/r + gradp[:,3]
    rotv = r*gradr[:,3] + 2*om - gradp[:,2]/r
    
    L = rho * np.power(s, 2.0/gam) / (2*rotv)

    print("    Plotting ")

    label = r"$\nabla \cdot v$"
    fig, ax = plt.subplots(1,1, figsize=(12,9))
    pde.primPlot(fig, ax, rjph, piph, r, divv, label, pars,
                    xlabel=r"$x$", ylabel=r"$y$",
                    planets=planetDat)
    plotname = "plot_eq_{0:s}_lin_div_v.png".format(name)
    print("    Saving "+plotname)
    fig.savefig(plotname)
    plt.close(fig)

    label = r"$\nabla \times v$"
    fig, ax = plt.subplots(1,1, figsize=(12,9))
    pde.primPlot(fig, ax, rjph, piph, r, rotv, label, pars,
                    xlabel=r"$x$", ylabel=r"$y$",
                    planets=planetDat)
    plotname = "plot_eq_{0:s}_lin_rot_v.png".format(name)
    print("    Saving "+plotname)
    fig.savefig(plotname)
    plt.close(fig)

    label = r"$L$"
    fig, ax = plt.subplots(1,1, figsize=(12,9))
    pde.primPlot(fig, ax, rjph, piph, r, L, label, pars,
                    xlabel=r"$x$", ylabel=r"$y$",
                    planets=planetDat)
    plotname = "plot_eq_{0:s}_lin_L.png".format(name)
    print("    Saving "+plotname)
    fig.savefig(plotname)
    plt.close(fig)

    label = r"$s$"
    fig, ax = plt.subplots(1,1, figsize=(12,9))
    pde.primPlot(fig, ax, rjph, piph, r, s, label, pars,
                    xlabel=r"$x$", ylabel=r"$y$",
                    planets=planetDat)
    plotname = "plot_eq_{0:s}_lin_s.png".format(name)
    print("    Saving "+plotname)
    fig.savefig(plotname)
    plt.close(fig)

    #for q in range(nq):
    #    label = r"$\partial_\phi$ " + vartex[q]
    #    fig, ax = plt.subplots(1,1, figsize=(12,9))
    #    pde.primPlot(fig, ax, rjph, piph, r, gradp[:,q], label, pars,
    #                    xlabel=r"$x$", ylabel=r"$y$",
    #                    planets=planetDat)
    #    plotname = "plot_eq_{0:s}_lin_gradp_{1:s}.png".format(name,varnames[q])
    #    print("    Saving "+plotname)
    #    fig.savefig(plotname)
    #    plt.close(fig)
    #    label = r"$\partial_r$ " + vartex[q]
    #    fig, ax = plt.subplots(1,1, figsize=(12,9))
    #    pde.primPlot(fig, ax, rjph, piph, r, gradr[:,q], label, pars,
    #                    xlabel=r"$x$", ylabel=r"$y$",
    #                    planets=planetDat)
    #    plotname = "plot_eq_{0:s}_lin_gradr_{1:s}.png".format(name,varnames[q])
    #    print("    Saving "+plotname)
    #    fig.savefig(plotname)
    #    plt.close(fig)

    fig, ax = plt.subplots(2,2)
    j = rjph.shape[0]/2
    ia = id0[j]
    ib = ia + nphi[j]
    f = prim[ia:ib,0]
    df = (f[2:] - f[:-2]) / (phi[ia:ib][2:] - phi[ia:ib][:-2])
    ax[0,0].plot(phi[ia:ib], f-f.mean(), 'k+')
    ax[0,0].plot(phi[ia:ib], gradp[ia:ib,0], 'b+')
    ax[0,0].plot(phi[ia:ib][1:-1], df, 'r+')
    f = prim[ia:ib,1]
    df = (f[2:] - f[:-2]) / (phi[ia:ib][2:] - phi[ia:ib][:-2])
    ax[0,1].plot(phi[ia:ib], f-f.mean(), 'k+')
    ax[0,1].plot(phi[ia:ib], gradp[ia:ib,1], 'b+')
    ax[0,1].plot(phi[ia:ib][1:-1], df, 'r+')
    f = prim[ia:ib,2]
    df = (f[2:] - f[:-2]) / (phi[ia:ib][2:] - phi[ia:ib][:-2])
    ax[1,0].plot(phi[ia:ib], f-f.mean(), 'k+')
    ax[1,0].plot(phi[ia:ib], gradp[ia:ib,2], 'b+')
    ax[1,0].plot(phi[ia:ib][1:-1], df, 'r+')
    f = prim[ia:ib,3]
    df = (f[2:] - f[:-2]) / (phi[ia:ib][2:] - phi[ia:ib][:-2])
    ax[1,1].plot(phi[ia:ib], f-f.mean(), 'k+')
    ax[1,1].plot(phi[ia:ib], gradp[ia:ib,3], 'b+')
    ax[1,1].plot(phi[ia:ib][1:-1], df, 'r+')
    plotname = "plot_phi_{0:s}.png".format(name)
    fig.savefig(plotname)
    plt.close(fig)
    
    fig, ax = plt.subplots(2,2)
    R = np.unique(r)[1:-1]
    f = prim[:,0]
    df = np.zeros(R.shape)
    ax[0,0].plot(r, f, 'k+')
    ax[0,0].plot(r, gradr[:,0], 'b+')
    ax[0,0].plot(R, df, 'r+')
    R = np.unique(r)[1:-1]
    f = prim[:,1]
    df = np.zeros(R.shape)
    ax[0,1].plot(r, f, 'k+')
    ax[0,1].plot(r, gradr[:,1], 'b+')
    ax[0,1].plot(R, df, 'r+')
    R = np.unique(r)[1:-1]
    f = prim[:,2]
    df = np.zeros(R.shape)
    ax[1,0].plot(r, f, 'k+')
    ax[1,0].plot(r, gradr[:,2], 'b+')
    ax[1,0].plot(R, df, 'r+')
    R = np.unique(r)[1:-1]
    f = prim[:,3]
    df = np.zeros(R.shape)
    ax[1,1].plot(r, f, 'k+')
    ax[1,1].plot(r, gradr[:,3], 'b+')
    ax[1,1].plot(R, df, 'r+')
    plotname = "plot_r_{0:s}.png".format(name)
    fig.savefig(plotname)
    plt.close(fig)

def calcCellGrad(grid, prim, PLM, PHI_MAX):

    rjph = grid[0]
    piph = grid[1]
    id0 = grid[2]
    nphi = grid[3]

    nr = rjph.shape[0]-1
    gradr = np.zeros(prim.shape, dtype=np.float)
    gradp = np.empty(prim.shape, dtype=np.float)
    dphi = np.empty(piph.shape, dtype=np.float)

    # Calculate dphi
    print("      Calc dphi")
    for j in xrange(nr):
        piphj = piph[id0[j]:id0[j]+nphi[j]]
        dphij = piphj - np.roll(piphj,1)
        while (dphij>PHI_MAX).any():
            dphij[dphij>PHI_MAX] -= PHI_MAX
        while (dphij<0.0).any():
            dphij[dphij<0.0] += PHI_MAX
        dphi[id0[j]:id0[j]+nphi[j]] = dphij[:]

    # Calculate phi gradients
    print("      Calc gp")
    for j in xrange(nr):
        ia = id0[j]
        ib = id0[j]+nphi[j]
        primC = prim[ia:ib,:]
        piphj = piph[ia:ib]
        dphij = dphi[ia:ib]
        
        phiC = piphj - 0.5*dphij
        phiR = np.roll(phiC,-1,axis=0)
        phiL = np.roll(phiC, 1,axis=0)
        primR = np.roll(primC,-1,axis=0)
        primL = np.roll(primC, 1,axis=0)

        #Make sure phiL,phiR on same branch as phiC
        while (phiL>phiC).any():
            phiL[phiL>phiC] -= PHI_MAX
        while (phiL<phiC-PHI_MAX).any():
            phiL[phiL<phiC-PHI_MAX] += PHI_MAX
        while (phiC>phiR).any():
            phiR[phiC>phiR] += PHI_MAX
        while (phiC<phiR-PHI_MAX).any():
            phiR[phiC<phiR-PHI_MAX] -= PHI_MAX

        gR = (primR[:,:] - primC[:,:]) / (phiR[:,None]-phiC[:,None])
        gC = (primR[:,:] - primL[:,:]) / (phiR[:,None]-phiL[:,None])
        gL = (primC[:,:] - primL[:,:]) / (phiC[:,None]-phiL[:,None])

        gradp[ia:ib,:] = minmod3(PLM*gL, gC, PLM*gR)

    dA = np.zeros(piph.shape)

    # Calc gradr - 1st pass
    for j in range(nr-1):
        j1 = j
        j2 = j+1
        i1a = id0[j1]
        i1b = i1a + nphi[j1]
        i2a = id0[j2]
        i2b = i2a + nphi[j2]

        rfm = rjph[j]
        rf = rjph[j+1]
        rfp = rjph[j+2]

        r1 = 2.0/3.0  * (rfm*rfm+rfm*rf+rf*rf) / (rfm+rf)
        r2 = 2.0/3.0  * (rf*rf+rf*rfp+rfp*rfp) / (rf+rfp)
        
        piph1 = piph[i1a:i1b]
        piph2 = piph[i2a:i2b]
        dphi1 = dphi[i1a:i1b]
        dphi2 = dphi[i2a:i2b]

        prim1 = prim[i1a:i1b,:]
        prim2 = prim[i2a:i2b,:]
        gradp1 = gradp[i1a:i1b,:]
        gradp2 = gradp[i2a:i2b,:]

        phi1p = piph1[0]
        phi1m = piph1[0] - dphi1[0]
        for i2 in range(nphi[j2]):
            phi2p = piph2[i2]
            while phi2p > phi1m+0.5*PHI_MAX:
                phi2p -= PHI_MAX
            while phi2p < phi1m-0.5*PHI_MAX:
                phi2p += PHI_MAX
            if phi2p > phi1m and phi2p - dphi2[i2] < phi1m:
                break

        for i1 in range(nphi[j1]):
            phi1p = piph1[i1]
            phi1m = piph1[i1] - dphi1[i1]
            phi1 = piph1[i1] - 0.5*dphi[i1]
            f1 = prim1[i1,:]
            while True:
                phi2p = piph2[i2]
                while phi2p > phi1m+0.5*PHI_MAX:
                    phi2p -= PHI_MAX
                while phi2p < phi1m-0.5*PHI_MAX:
                    phi2p += PHI_MAX
                phi2m = phi2p - dphi2[i2]
                phi2 = phi2p - 0.5*dphi2[i2]

                f2 = prim2[i2,:]

                phip = min(phi1p, phi2p)
                phim = max(phi1m, phi1m)

                dAf = rf*(phip-phim)
                phif = 0.5*(phim+phip)

                s = (f2 + gradp2[i2]*(phif-phi2) - f1 - gradp1[i1]*(phif-phi1)
                        ) / (r2 - r1)

                gradr[i1a+i1,:] += s[:]*dAf
                gradr[i2a+i2,:] += s[:]*dAf
                dA[i1a+i1] += dAf
                dA[i2a+i2] += dAf

                if phi2p > phi1p:
                    break
                else:
                    i2 += 1
                    if i2 >= nphi[j2]:
                        i2 = 0;

    gradr[:,:] /= dA[:,None]

    # Calc gradr - 2nd pass
    for j in range(nr-1):
        j1 = j
        j2 = j+1
        i1a = id0[j1]
        i1b = i1a + nphi[j1]
        i2a = id0[j2]
        i2b = i2a + nphi[j2]

        rfm = rjph[j]
        rf = rjph[j+1]
        rfp = rjph[j+2]

        r1 = 2.0/3.0  * (rfm*rfm+rfm*rf+rf*rf) / (rfm+rf)
        r2 = 2.0/3.0  * (rf*rf+rf*rfp+rfp*rfp) / (rf+rfp)

        piph1 = piph[i1a:i1b]
        piph2 = piph[i2a:i2b]
        dphi1 = dphi[i1a:i1b]
        dphi2 = dphi[i2a:i2b]

        prim1 = prim[i1a:i1b,:]
        prim2 = prim[i2a:i2b,:]
        gradp1 = gradp[i1a:i1b,:]
        gradp2 = gradp[i2a:i2b,:]

        phi1p = piph1[0]
        phi1m = piph1[0] - dphi1[0]
        for i2 in range(nphi[j2]):
            phi2p = piph2[i2]
            while phi2p > phi1m+0.5*PHI_MAX:
                phi2p -= PHI_MAX
            while phi2p < phi1m-0.5*PHI_MAX:
                phi2p += PHI_MAX
            if phi2p > phi1m and phi2p - dphi2[i2] < phi1m:
                break

        for i1 in range(nphi[j1]):
            phi1p = piph1[i1]
            phi1m = piph1[i1] - dphi1[i1]
            phi1 = piph1[i1] - 0.5*dphi[i1]
            f1 = prim1[i1,:]
            while True:
                phi2p = piph2[i2]
                while phi2p > phi1m+0.5*PHI_MAX:
                    phi2p -= PHI_MAX
                while phi2p < phi1m-0.5*PHI_MAX:
                    phi2p += PHI_MAX
                phi2m = phi2p - dphi2[i2]
                phi2 = phi2p - 0.5*dphi2[i2]

                f2 = prim2[i2,:]

                phip = min(phi1p, phi2p)
                phim = max(phi1m, phi1m)

                dAf = rf*(phip-phim)
                phif = 0.5*(phim+phip)

                s = (f2 + gradp2[i2]*(phif-phi2) - f1 - gradp1[i1]*(phif-phi1)
                        ) / (r2 - r1)

                gradr[i1a+i1,:] = minmod2(gradr[i1a+i1,:], PLM*s)
                gradr[i2a+i2,:] = minmod2(gradr[i2a+i2,:], PLM*s)

                if phi2p > phi1p:
                    break
                else:
                    i2 += 1
                    if i2 >= nphi[j2]:
                        i2 = 0;

    return gradr, gradp

def minmod2(a, b):
    res = np.zeros(a.shape)
    altb = np.fabs(a) <= np.fabs(b)
    isame = (a*b > 0)
    res[altb*isame] = a[altb*isame]
    res[(-altb)*isame] = b[(-altb)*isame]

    return res

def minmod3(a, b, c):
    res = np.zeros(a.shape)
    altb = np.fabs(a) <= np.fabs(b)
    bltc = np.fabs(b) <= np.fabs(c)
    clta = np.fabs(c) <= np.fabs(a)
    isame = ((a*b)>0) * ((b*c)>0)
    ia = altb * (-clta) * isame
    ib = bltc * (-altb) * isame
    ic = clta * (-bltc) * isame

    res[ia] = a[ia]
    res[ib] = b[ib]
    res[ic] = c[ic]

    return res



if __name__ == "__main__":

    parser = ap.ArgumentParser(description="Create 2D plots of gradients.")
    parser.add_argument('checkpoints', nargs='+', 
                            help="Checkpoint (.h5) files to plot.")

    args = parser.parse_args()

    files = args.checkpoints

    for f in files:
        plotGradients(f)

