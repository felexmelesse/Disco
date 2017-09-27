import sys
import math
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import discoUtil as du

xscale = "log"
yscale = "log"
GAM = 1.66666666667
RMAX = 0.8
GR = False

def plotCheckpoint(file):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = du.loadCheckpoint(file)

    rho = prim[:,0]
    P = prim[:,1]
    ur = prim[:,2]
    up = prim[:,3]
    uz = prim[:,4]
    Br = prim[:,5]
    Bp = prim[:,6]
    Bz = prim[:,7]

    piph = dat[3]
    rjph = dat[0]
    zkph = dat[1]

    B2 = Br*Br + Bp*Bp + Bz*Bz
    u2 = ur*ur + r*r*up*up + uz*uz
    w2 = 1.0 + u2
    uB = ur*Br + r*up*Bp + uz*Bz
    b2 = (B2 + uB*uB) / w2
    rhoh = rho + GAM/(GAM-1.0)*P

    cs = np.sqrt(GAM*P/rhoh)
    cA = np.sqrt(b2/(rhoh+b2))
    Ma = np.sqrt(u2) / (cs/np.sqrt(1-cs*cs))

    NR = len(rjph)-1
    NZ = len(zkph)-1

    rF = np.empty((NZ, NR))
    zF = np.empty((NZ, NR))
    drF = np.empty((NZ, NR))
    dzF = np.empty((NZ, NR))
    rhoF = np.empty((NZ, NR))
    PF = np.empty((NZ, NR))
    urF = np.empty((NZ, NR))
    upF = np.empty((NZ, NR))
    uzF = np.empty((NZ, NR))
    BrF = np.empty((NZ, NR))
    BpF = np.empty((NZ, NR))
    BzF = np.empty((NZ, NR))
    
    RrpF = np.empty((NZ, NR))
    MrpF = np.empty((NZ, NR))
    PbF  = np.empty((NZ, NR))

    vol = 0.0
    for k in range(NZ):
        dz = zkph[k+1] - zkph[k]
        indz = (z > zkph[k]) * (z < zkph[k+1])
        for j in range(NR):
            dr = rjph[j+1] - rjph[j]
            indr = (r > rjph[j]) * (r < rjph[j+1])
            ind = indr * indz
            piphi = piph[ind]
            dphi = piphi - np.roll(piphi, 1)
            dphi[dphi < 0] = dphi[dphi>0].mean()

            DP = dphi.sum()

            rF[k,j]   = 0.5*(r[j] + r[j+1])
            zF[k,j]   = 0.5*(z[k] + z[k+1])
            drF[k,j]  = dr
            dzF[k,j]  = dz
            rhoF[k,j] = (rho[ind]*dphi).mean() / DP
            PF[k,j]   =   (P[ind]*dphi).mean() / DP
            urF[k,j]  =  (ur[ind]*dphi).mean() / DP
            upF[k,j]  =  (up[ind]*dphi).mean() / DP
            uzF[k,j]  =  (uz[ind]*dphi).mean() / DP
            BrF[k,j]  =  (Br[ind]*dphi).mean() / DP
            BpF[k,j]  =  (Bp[ind]*dphi).mean() / DP
            BzF[k,j]  =  (Bz[ind]*dphi).mean() / DP
            
            if not GR:
                RrpF[k,j] =  (rho[ind]*ur[ind]*rF[k,j]*(up[ind]-upF[k,j])*dphi
                            ).mean() / DP
            else:
                RrpF[k,j] =  (rho[ind]*ur[ind]*(up[ind]-upF[k,j])/rF[k,j]*dphi
                            ).mean() / DP
            MrpF[k,j] =  (-Br[ind]*Bp[ind]*dphi).mean() / DP
            PbF[k,j]   =  (0.5*B2[ind]*dphi).mean() / DP

            vol += rF[k,j] * dr * dz * DP
            

    Rrp = (rF * RrpF * drF * dzF).sum() / vol
    Mrp = (rF * MrpF * drF * dzF).sum() / vol
    Pg  = (rF * PF   * drF * dzF).sum() / vol
    Pb  = (rF * PbF  * drF * dzF).sum() / vol

    print("al:   {0:.6g}".format((Rrp + Mrp) / Pg))
    print("al_M: {0:.6g}".format(Mrp / Pg))
    print("be:   {0:.6g}".format(Pg / Pb))
    print("thB:  {0:.6g}".format(180.0/math.pi * 0.5*math.asin(Mrp/Pb)))

    print("   Plotting...")
    nq = prim.shape[1]

    title = "DISCO t = {0:.3g}".format(t)

    fig, ax = plt.subplots(2,2)
    C = ax[0,0].pcolormesh(rjph, zkph, RrpF, cmap=plt.cm.inferno)
    fig.colorbar(C, ax=ax[0,0])
    ax[0,0].set_title(r"$R_{r\phi}$")
    ax[0,0].set_ylabel(r"$z$")
    C = ax[0,1].pcolormesh(rjph, zkph, MrpF, cmap=plt.cm.inferno)
    fig.colorbar(C, ax=ax[0,1])
    ax[0,1].set_title(r"$M_{r\phi}$")
    C = ax[1,0].pcolormesh(rjph, zkph, PF, cmap=plt.cm.inferno)
    fig.colorbar(C, ax=ax[1,0])
    ax[1,0].set_title(r"$P_g$")
    ax[1,0].set_xlabel(r"$r$")
    ax[1,0].set_ylabel(r"$z$")
    C = ax[1,1].pcolormesh(rjph, zkph, PbF, cmap=plt.cm.inferno)
    fig.colorbar(C, ax=ax[1,1])
    ax[1,1].set_title(r"$P_b$")
    ax[1,1].set_xlabel(r"$r$")

    fig.suptitle(title, fontsize=18)

    plt.tight_layout()

    name = file.split('/')[-1].split('.')[0].split('_')[-1]
    plotname = "plot_mri_{0:s}.png".format(name)
    
    print("   Saving {0:s}...".format(plotname))
    fig.savefig(plotname)

    plt.close(fig)

    return t, (Rrp, Mrp, Pg, Pb)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Makes plots of Disco prims as a function of r.")
        print("usage: python plotDiscoR.py <checkpoint.h5 ...>")
        sys.exit()

    files = sys.argv[1:]
    dat = np.empty((len(files), 4))
    t = np.empty(len(files))
    for i,f in enumerate(files):
        ret = plotCheckpoint(f)
        t[i] = ret[0]
        dat[i,:] = np.array(ret[1])

    Rrp = dat[:,0]
    Mrp = dat[:,1]
    Pg = dat[:,2]
    Pb = dat[:,3]

    fig, ax = plt.subplots(2,3)
    ax[0,0].plot(t, Rrp, label=r"$R_{r\phi}$")
    ax[0,0].plot(t, Mrp, label=r"$M_{r\phi}$")
    ax[0,0].legend()
    ax[0,1].plot(t, Pg, label=r"$P_g$")
    ax[0,1].plot(t, Pb, label=r"$P_b$")
    ax[0,1].legend()
    ax[1,0].plot(t, (Rrp+Mrp)/Pg, label=r"$\alpha$")
    ax[1,0].plot(t, Mrp/Pg, label=r"$\alpha_M$")
    ax[1,0].legend()
    ax[1,1].plot(t, Pg/Pb, label=r"$\beta$")
    ax[1,1].legend()
    ax[1,2].plot(t, 180.0/np.pi * 0.5*np.arcsin(Mrp / Pb), label=r"$\theta_B$")
    ax[1,2].legend()
    fig.tight_layout()

    name = "mri"
    fig.savefig(name+".pdf")

