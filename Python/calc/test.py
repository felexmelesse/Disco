import numpy as np
import matplotlib.pyplot as plt
import _calc as ca

plotNewt = False
plotRel = True

xscale = 'log'
yscale = 'log'

gam = 4.0/3.0
Mdot = 1.0
rho0 = 1.0e-2
GM = 1.0
a0 = 1.7e-1
x0 = 3.0
x1 = 4.0e1

r = np.logspace(np.log10(x0*GM), np.log10(x1*GM), 200)

rhoN, uN, PN = ca.bondi_newt(Mdot, GM, gam, rho0, r)
aN = np.sqrt(gam*PN/rhoN)

rhoR, uR, PR = ca.bondi_rel(Mdot, GM, gam, a0, r)
aR = np.sqrt(gam*PR/(rhoR + gam*PR/(gam-1.0)))
u2R = uR*uR / (1.0-2*GM/r)
machR = np.sqrt(u2R*(1-aR*aR)) / aR

if plotNewt:

    fig1, ax1 = plt.subplots(2,2)

    ax1[0,0].plot(r, rhoN, 'k+')
    ax1[0,0].set_xscale(xscale)
    ax1[0,0].set_yscale(yscale)
    ax1[0,1].plot(r, PN, 'k+')
    ax1[0,1].set_xscale(xscale)
    ax1[0,1].set_yscale(yscale)
    ax1[1,0].plot(r, -uN, 'k+')
    ax1[1,0].set_xscale(xscale)
    ax1[1,0].set_yscale(yscale)
    ax1[1,1].plot(r, aN, 'k+')
    ax1[1,1].set_xscale(xscale)
    ax1[1,1].set_yscale(yscale)

if plotRel:

    fig2, ax2 = plt.subplots(3,3, figsize=(14,10))

    ax2[0,0].plot(r, rhoR, 'k+')
    ax2[0,0].set_ylabel(r"$\rho$")
    ax2[0,0].set_xscale(xscale)
    ax2[0,0].set_yscale(yscale)

    ax2[0,1].plot(r, PR, 'k+')
    ax2[0,1].set_ylabel(r"$P$")
    ax2[0,1].set_xscale(xscale)
    ax2[0,1].set_yscale(yscale)

    ax2[1,0].plot(r, -uR, 'k+')
    ax2[1,0].set_ylabel(r"$u^r$")
    ax2[1,0].set_xscale(xscale)
    ax2[1,0].set_yscale(yscale)

    ax2[1,1].plot(r, aR, 'k+')
    ax2[1,1].set_ylabel(r"$a$")
    ax2[1,1].set_xscale(xscale)
    ax2[1,1].set_yscale(yscale)

    ax2[1,2].plot(r, machR, 'k+')
    ax2[1,2].set_ylabel(r"$\mathcal{M}$")
    ax2[1,2].set_xscale(xscale)
    ax2[1,2].set_yscale(yscale)

    ax2[2,0].plot(r, -4*np.pi*r*r*rhoR*uR, 'k+')
    ax2[2,0].plot(r, Mdot*np.ones(r.shape), 'b')
    ax2[2,0].set_ylabel(r"$\dot{M}$")
    ax2[2,0].set_xscale(xscale)
    ax2[2,0].set_yscale("linear")

if plotNewt or plotRel:
    plt.show()


