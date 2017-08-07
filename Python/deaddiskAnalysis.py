import sys
import math
import argparse as ap
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discoUtil as du

def moment_arm(r1, r2):
    return math.sqrt(0.5*(r1*r1+r2*r2))

def calc_grav(r, phi, planetDat):

    gr = np.zeros(phi.shape)
    gp = np.zeros(phi.shape)

    x = r*np.cos(phi)
    y = r*np.sin(phi)

    npl = planetDat.shape[0]

    for p in range(npl):
        Mp = planetDat[p,0]
        rp = planetDat[p,3]
        pp = planetDat[p,4]
        xp = rp * math.cos(pp)
        yp = rp * math.sin(pp)
        dx = x-xp
        dy = y-yp
        rp = np.sqrt(dx*dx + dy*dy)

        gx = -Mp*dx/(rp*rp*rp)
        gy = -Mp*dy/(rp*rp*rp)

        gr[:] +=  (x*gx + y*gy)/r
        gp[:] += -y*gx + x*gy

    return gr, gp

def fftPhi(rpjh, prim, PHIMAX):
    return

def fluxPlots(grid, prim, planetDat, pars, name):

    rjph = grid[0]
    piph = grid[1]
    id0 = grid[2]
    nphi = grid[3]

    PHIMAX = pars['Phi_Max']

    nr = rjph.shape[0]-1
    N = prim.shape[0]
    nq = prim.shape[1]

    R = np.zeros((nr,))
    Sig = np.zeros((nr,))
    Mdot = np.zeros((nr,))
    L = np.zeros((nr,))
    FJ = np.zeros((nr,))
    Tg = np.zeros((nr,))

    for j in range(nr):
        ia = id0[j]
        ib = ia + nphi[j]

        r = moment_arm(rjph[j], rjph[j+1])
        dphi = piph[ia:ib] - np.roll(piph[ia:ib],1)
        while (dphi > PHIMAX).any():
            dphi[dphi>PHIMAX] -= PHIMAX
        while (dphi < 0.0).any():
            dphi[dphi<0.0] += PHIMAX
        phi = piph[ia:ib] - 0.5*dphi

        sig = prim[ia:ib,0]
        P = prim[ia:ib,1]
        vr = prim[ia:ib,2]
        om = prim[ia:ib,3]

        gr, gp = calc_grav(r, phi, planetDat)

        dA = r*dphi
        
        R[j] = r
        Sig[j] = (sig * dA).sum()
        Mdot[j] = -(sig*vr * dA).sum()
        L[j] = (sig*om*r*r * dA).sum()
        FJ[j] = (sig*vr*r*r*om * dA).sum()
        Tg[j] = (sig*gp * dA).sum()

    A = R*PHIMAX
    ell = L / Sig
    FRe = FJ + Mdot*ell

    TJ = np.zeros((nr,))
    Tacc = np.zeros((nr,))
    TRe = np.zeros((nr,))

    TJ[1:-1] = - (FJ[2:]-FJ[:-2]) / (R[2:]-R[:-2])
    Tacc[1:-1] = Mdot[1:-1] * (ell[2:]-ell[:-2])/(R[2:]-R[:-2])
    TRe[1:-1] = - (FRe[2:]-FRe[:-2]) / (R[2:]-R[:-2])

    fig, ax = plt.subplots(2,3, figsize=(16,9))
    ax[0,0].plot(R[2:-2], (Sig/A)[2:-2], 'k+')
    ax[0,0].set_xlabel(r'$r$')
    ax[0,0].set_ylabel(r'$\hat{\Sigma} = \langle \Sigma \rangle / 2\pi r$')
    ax[0,1].plot(R[2:-2], Mdot[2:-2], 'k+')
    ax[0,1].set_xlabel(r'$r$')
    ax[0,1].set_ylabel(r'$\dot{M} = -\langle \Sigma v_r \rangle$')
    ax[0,2].plot(R[2:-2], ell[2:-2], 'k+')
    ax[0,2].set_xlabel(r'$r$')
    ax[0,2].set_ylabel(
      r'$\ell = -\langle \Sigma r^2 \Omega \rangle / \langle \Sigma \rangle$')
    ax[1,0].plot(R[2:-2], FJ[2:-2], 'k+', 
                    label=r'$\langle \Sigma r^2 \Omega v_r \rangle$')
    ax[1,0].plot(R[2:-2], -(Mdot*ell)[2:-2], 'b+', label=r'$-\dot{M} \ell$')
    ax[1,0].plot(R[2:-2], FRe[2:-2], 'r+',
            label=r'$\langle \Sigma r^2 \Omega v_r \rangle_{\mathrm{Re}}$')
    ax[1,0].set_xlabel(r'$r$')
    ax[1,0].set_ylabel(r'$F_J$')
    ax[1,0].legend(loc='lower left', fontsize=12)
    
    ax[1,1].plot(R[2:-2], TJ[2:-2], 'k+', label=r'$\tau_{F_J}$')
    ax[1,1].plot(R[2:-2], Tg[2:-2], 'g+', label=r'$\tau_g$')
    ax[1,1].set_xlabel(r'$r$')
    ax[1,1].set_ylabel(r'$\tau$')
    ax[1,1].legend(loc='lower left', fontsize=12)

    ax[1,2].plot(R[2:-2], Tacc[2:-2], 'b+', label=r'$\tau_{\dot{M}}$')
    ax[1,2].plot(R[2:-2], TRe[2:-2], 'r+', label=r'$\tau_{\mathrm{Re}}$')
    ax[1,2].plot(R[2:-2], Tg[2:-2], 'g+', label=r'$\tau_g$')
    ax[1,2].set_xlabel(r'$r$')
    ax[1,2].set_ylabel(r'$\tau$')
    ax[1,2].legend(loc='lower left', fontsize=12)

    fig.tight_layout()

    plotname = "plot_dd_flux_{0:s}.pdf".format(name)
    fig.savefig(plotname)
    plt.close(fig)

    return R, Sig, ell, Mdot, FJ, Tg, TJ, Tacc, TRe

def fourierPlots(grid, prim, planetDat, pars, name, plot_heatmap=False,
                    nm=32):

    rjph = grid[0]
    piph = grid[1]
    id0 = grid[2]
    nphi = grid[3]

    PHIMAX = pars['Phi_Max']
    gam = pars['Adiabatic_Index']
    R0 = pars['Init_Par1']
    DR = pars['Init_Par2']
    DR2 = pars['Init_Par3']
    R1 = R0-DR-DR2/2
    R2 = R0+DR+DR2/2

    nr = rjph.shape[0]-1
    nq = prim.shape[1]

    R = np.zeros((nr,))

    if nm is None:
        nm = min([n/2+1 for n in nphi])

    sigm = np.empty((nm,nr), dtype=np.complex)
    Pm = np.empty((nm,nr), dtype=np.complex)
    vrm = np.empty((nm,nr), dtype=np.complex)
    omm = np.empty((nm,nr), dtype=np.complex)
    hm = np.empty((nm,nr), dtype=np.complex)
    grm = np.empty((nm,nr), dtype=np.complex)
    gpm = np.empty((nm,nr), dtype=np.complex)


    for j in range(nr):

        ia = id0[j]
        ib = ia + nphi[j]

        r = moment_arm(rjph[j], rjph[j+1])

        dphi = piph[ia:ib] - np.roll(piph[ia:ib],1)
        while (dphi > PHIMAX).any():
            dphi[dphi>PHIMAX] -= PHIMAX
        while (dphi < 0.0).any():
            dphi[dphi<0.0] += PHIMAX
        phi = piph[ia:ib] - 0.5*dphi

        gr, gp = calc_grav(r, phi, planetDat)

        sig = prim[ia:ib,0]
        P = prim[ia:ib,1]
        vr = prim[ia:ib,2]
        om = prim[ia:ib,3]

        h = gam/(gam-1) * P/sig

        R[j] = r
        sigm[:,j] = np.fft.rfft(sig)[:nm]
        Pm[:,j] = np.fft.rfft(P)[:nm]
        vrm[:,j] = np.fft.rfft(vr)[:nm]
        omm[:,j] = np.fft.rfft(om)[:nm]
        hm[:,j] = np.fft.rfft(h)[:nm]
        grm[:,j] = np.fft.rfft(gr)[:nm]
        gpm[:,j] = np.fft.rfft(gp)[:nm]

    Mb = np.arange(nm+1) - 0.5

    if plot_heatmap:
        fig, ax = plt.subplots(2, 4, figsize=(16,9))
        C = ax[0,0].pcolormesh(rjph[2:-2], Mb, np.abs(sigm)[:,2:-2], 
                            cmap=mpl.cm.viridis)
        ax[0,0].set_ylabel(r'$m$')
        ax[0,0].set_title(r'$\Sigma$')
        fig.colorbar(C, ax=ax[0,0])
        

        C = ax[0,1].pcolormesh(rjph[2:-2], Mb, np.abs(Pm)[:,2:-2], 
                cmap=mpl.cm.viridis)
        ax[0,1].set_ylabel(r'$m$')
        ax[0,1].set_title(r'$P$')
        fig.colorbar(C, ax=ax[0,1])

        C = ax[0,2].pcolormesh(rjph[2:-2], Mb, np.abs(vrm)[:,2:-2], 
                cmap=mpl.cm.viridis)
        ax[0,2].set_ylabel(r'$m$')
        ax[0,2].set_title(r'$v_r$')
        fig.colorbar(C, ax=ax[0,2])

        C = ax[0,3].pcolormesh(rjph[2:-2], Mb, np.abs(omm)[:,2:-2], 
                cmap=mpl.cm.viridis)
        ax[0,3].set_ylabel(r'$m$')
        ax[0,3].set_title(r'$\Omega$')
        fig.colorbar(C, ax=ax[0,3])

        C = ax[1,0].pcolormesh(rjph[2:-2], Mb, np.abs(hm)[:,2:-2], 
                            cmap=mpl.cm.viridis)
        ax[1,0].set_xlabel(r'$r$')
        ax[1,0].set_ylabel(r'$m$')
        ax[1,0].set_title(r'$h$')
        fig.colorbar(C, ax=ax[1,0])

        C = ax[1,1].pcolormesh(rjph[2:-2], Mb, np.abs(grm)[:,2:-2], 
                cmap=mpl.cm.viridis)
        ax[1,1].set_xlabel(r'$r$')
        ax[1,1].set_ylabel(r'$m$')
        ax[1,1].set_title(r'$g_r$')
        fig.colorbar(C, ax=ax[1,1])

        C = ax[1,2].pcolormesh(rjph[2:-2], Mb, np.abs(gpm)[:,2:-2], 
                cmap=mpl.cm.viridis)
        ax[1,2].set_xlabel(r'$r$')
        ax[1,2].set_ylabel(r'$m$')
        ax[1,2].set_title(r'$g_\phi$')
        fig.colorbar(C, ax=ax[1,2])

        fig.tight_layout()

        plotname = "plot_dd_fourier_{0:s}.png".format(name)
        fig.savefig(plotname)
        plt.close(fig)

    fig, ax = plt.subplots(2,4, figsize=(16,9))
    for m in range(1,9):
        ax[0,0].plot(R, sigm[m,:].real, label=r'$m={0}$'.format(m))
        ax[0,1].plot(R, Pm[m,:].real)
        ax[0,2].plot(R, vrm[m,:].real)
        ax[0,3].plot(R, omm[m,:].real)
        ax[1,0].plot(R, hm[m,:].real)
        ax[1,1].plot(R, grm[m,:].real)
        ax[1,2].plot(R, gpm[m,:].real)
    ax[0,0].legend(fontsize=10)
    ax[0,0].set_ylabel(r"$\Sigma$")
    ax[0,1].set_ylabel(r"$P$")
    ax[0,2].set_ylabel(r"$v_r$")
    ax[0,3].set_ylabel(r"$\Omega$")
    ax[1,0].set_ylabel(r"$h$")
    ax[1,1].set_ylabel(r"$g_r$")
    ax[1,2].set_ylabel(r"$g_\phi$")
    fig.tight_layout()

    fig, ax = plt.subplots(2,4, figsize=(16,9))
    ind = (R>R1)*(R<R2)
    for m in range(1,9):
        ax[0,0].plot(R[ind], sigm[m,ind].real, label=r'$m={0}$'.format(m))
        ax[0,1].plot(R[ind], Pm[m,ind].real)
        ax[0,2].plot(R[ind], vrm[m,ind].real)
        ax[0,3].plot(R[ind], omm[m,ind].real)
        ax[1,0].plot(R[ind], hm[m,ind].real)
        ax[1,1].plot(R[ind], grm[m,ind].real)
        ax[1,2].plot(R[ind], gpm[m,ind].real)
    ax[0,0].legend(fontsize=10)
    ax[0,0].set_ylabel(r"$\Sigma$")
    ax[0,1].set_ylabel(r"$P$")
    ax[0,2].set_ylabel(r"$v_r$")
    ax[0,3].set_ylabel(r"$\Omega$")
    ax[1,0].set_ylabel(r"$h$")
    ax[1,1].set_ylabel(r"$g_r$")
    ax[1,2].set_ylabel(r"$g_\phi$")
    fig.tight_layout()

    plotname = "plot_dd_fourier_real_cut_{0:s}.pdf".format(name)
    fig.savefig(plotname)
    plt.close(fig)

    plotname = "plot_dd_fourier_real_{0:s}.pdf".format(name)
    fig.savefig(plotname)
    plt.close(fig)
    
    fig, ax = plt.subplots(2,4, figsize=(16,9))
    ax[0,0].plot(R, sigm[0,:].real)
    ax[0,1].plot(R, Pm[0,:].real)
    ax[0,2].plot(R, vrm[0,:].real)
    ax[0,3].plot(R, omm[0,:].real)
    ax[1,0].plot(R, hm[0,:].real)
    ax[1,1].plot(R, grm[0,:].real)
    ax[1,2].plot(R, gpm[0,:].real)
    ax[0,0].set_ylabel(r"$\Sigma$")
    ax[0,1].set_ylabel(r"$P$")
    ax[0,2].set_ylabel(r"$v_r$")
    ax[0,3].set_ylabel(r"$\Omega$")
    ax[1,0].set_ylabel(r"$h$")
    ax[1,1].set_ylabel(r"$g_r$")
    ax[1,2].set_ylabel(r"$g_\phi$")
    fig.tight_layout()

    plotname = "plot_dd_aveR_{0:s}.pdf".format(name)
    fig.savefig(plotname)
    plt.close(fig)

    return sigm, Pm, vrm, omm, hm, grm, gpm

def summaryQuantities(grid, prim, pars, planetDat, name):

    rjph = grid[0]
    piph = grid[1]
    id0 = grid[2]
    nphi = grid[3]

    PHIMAX = pars['Phi_Max']
    gam = pars['Adiabatic_Index']
    R0 = pars['Init_Par1']
    DR = pars['Init_Par2']
    DR2 = pars['Init_Par3']
    R1 = R0-DR-DR2/2
    R2 = R0+DR+DR2/2

    nr = rjph.shape[0]-1
    nq = prim.shape[1]
    N = prim.shape[0]

    R = np.zeros((nr,))

    M = 0.0
    Mq = 0.0

    dphi = np.empty((N,))
    r = np.empty((N,))

    for j in range(nr):
        ia = id0[j]
        ib = ia + nphi[j]

        rj = moment_arm(rjph[j], rjph[j+1])
        dphij = piph[ia:ib] - np.roll(piph[ia:ib],1)
        while (dphij > PHIMAX).any():
            dphij[dphij>PHIMAX] -= PHIMAX
        while (dphij < 0.0).any():
            dphij[dphij<0.0] += PHIMAX
        
        r[ia:ib] = rj
        dphi[ia:ib] = dphij[:]

    dV = r*dphi

    sig = prim[:,0]
    P = prim[:,1]
    vr = prim[:,2]
    om = prim[:,3]
    q = prim[:,5]
    s = 1.0/(gam-1.0) * np.log(P * np.power(sig,-gam))

    gr, gp = calc_grav(r, phi, planetDat)

    M = (sig * dV).sum()
    Mq = (q * sig * dV).sum()
    L = (sig*r*r*om * dV).sum()
    Lq = (q * sig*r*r*om * dV).sum()
    Eth = (P/(gam-1) * dV).sum()
    Ethq = (q * P/(gam-1) * dV).sum()
    Ek = (0.5*sig*(vr*vr+r*r*om*om) * dV).sum()
    Ekq = (q * 0.5*sig*(vr*vr+r*r*om*om) * dV).sum()
    S = (sig*s * dV).sum()
    Sq = (a * sig*s * dV).sum()
    

    return M, Mq, L, Lq, Eth, Ethq, Ek, Ekq, S, Sq

def analysisSingle(filename):

    print("Loading {0}".format(filename))
    t, r, phi, z, prim, dat = du.loadCheckpoint(filename)
    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    piph = dat[3]
    planetDat = dat[4]
    nphi = dat[7]
    pars = du.loadPars(filename)

    k = int(zkph.shape[0]/2-1)
    ind = (z>zkph[k]) * (z<zkph[k+1])

    r = r[ind]
    phi = phi[ind]
    z = z[ind]
    prim = prim[ind,:]
    piph = piph[ind]
    primPhi0 = primPhi0[k,:]
    nphi = nphi[k,:]
    id0 = np.zeros(nphi.shape, dtype=np.int)
    id0[1:] = np.cumsum(nphi)[:-1]

    grid = (rjph, piph, id0, nphi)

    name = filename.split('/')[-1].split('.')[0].split('_')[-1]

    print("    Fluxes...")
    fluxDat = fluxPlots(grid, prim, planetDat, pars, name)
    R, Sig, ell, Mdot, FJ, Tg, TJ, Tacc, TRe = fluxDat

    print("    Fourier...")
    fourierPlots(grid, prim, planetDat, pars, name)
    print("    Summary...")
    summ = summaryQuantities(grid, prim, pars, planetDat, name)

    return t, summ


def analysisSumm(summaries):

    N = len(summaries)

    T = np.zeros(N)
    M = np.zeros(N)
    Mq = np.zeros(N)
    L = np.zeros(N)
    Lq = np.zeros(N)
    Eth = np.zeros(N)
    Ethq = np.zeros(N)
    Ek = np.zeros(N)
    Ekq = np.zeros(N)
    S = np.zeros(N)
    Sq = np.zeros(N)

    for i,summ in enumerate(summaries):
        t = summ[0]
        dat = summ[1]
        T[i] = t
        M[i], Mq[i], L[i], Lq[i], Eth[i], Ethq[i], Ek[i], Ekq[i], S[i], Sq[i] = dat

    fig, ax = plt.subplots(2,2, figsize=(12,9))
    ax[0,0].plot(T, M, 'k+')
    ax[0,0].plot(T, Mq, 'b+')
    ax[0,0].set_ylabel(r'$M$')
    
    ax[0,1].plot(T, L, 'k+')
    ax[0,1].plot(T, Lq, 'b+')
    ax[0,1].set_ylabel(r'$J$')
    
    ax[1,0].plot(T, Eth, 'k+')
    ax[1,0].plot(T, Ethq, 'b+')
    ax[1,0].set_xlabel(r'$t$')
    ax[1,0].set_ylabel(r'$J$')
    
    ax[1,1].plot(T, S, 'k+')
    ax[1,1].plot(T, Sq, 'b+')
    ax[1,1].set_xlabel(r'$t$')
    ax[1,1].set_ylabel(r'$J$')

    plotname = "plot_dd_summ_cons_time.pdf"
    fig.savefig(plotname)
    plt.close(fig)


if __name__ == "__main__":

    parser = ap.ArgumentParser(description="Perform deaddisk analysis")
    parser.add_argument('checkpoints', nargs='+', 
                            help="Checkpoint (.h5) files to plot.")

    args = parser.parse_args()
    files = args.checkpoints

    summaries = []
    for f in files:
        summ = analysisSingle(f)
        summaries.append(summ)

    analysisSumm(summaries)
