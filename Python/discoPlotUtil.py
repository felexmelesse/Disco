import os
import sys
import argparse as ag
import h5py as h5
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discoUtil as du
import discoGeom as dg

def plotZSlice(fig, ax, rjph, piph, r, q, Z, label, pars, opts, vmin=None, 
            vmax=None, noGhost=False, colorbar=True, xlabel=None, ylabel=None,
            log=False, rmax=None, planets=None, cmap=None):

    phi_max = pars['Phi_Max']

    if cmap is None:
        try:
            cmap = mpl.cm.inferno
        except:
            cmap = mpl.cm.afmhot
        
    if vmin is None:
        vmin = q.min()
    if vmax is None:
        vmax = q.max()

    Rs = np.unique(r)
    if noGhost:
        Rs = Rs[:-2]

    if rmax is None or rmax <=0.0:
        lim_float = True
    else:
        lim_float = False

    if log:
        norm = mpl.colors.LogNorm(vmin, vmax)
    else:
        norm = mpl.colors.Normalize(vmin, vmax)

    pimh_min = np.inf

    for i, R in enumerate(Rs):
        ind = r==R
        imax = np.argmax(piph[ind])

        apiph = np.roll(piph[ind], -imax-1)
        aq = np.roll(q[ind], -imax-1)

        phif = np.empty(len(apiph)+1)
        phif[1:] = apiph[:]
        phif[0] = apiph[-1] - phi_max

        rf = rjph[i:i+2]

        #x = rf[:,None] * np.cos(phif)[None,:]
        #y = rf[:,None] * np.sin(phif)[None,:]

        x, y, z = dg.getXYZ(rf[:,None], phif[None,:], Z.mean(), opts, pars)

        C = ax.pcolormesh(x, y, aq[None,:], 
                cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)

        if lim_float and rf.max() > rmax:
            rmax = rf.max()
        if phif.min() < pimh_min:
            pimh_min =  phif.min()

    if planets is not None:
        rpl = planets[:,3]
        ppl = planets[:,4]
        xpl = rpl * np.cos(ppl)
        ypl = rpl * np.sin(ppl)
        ax.plot(xpl, ypl, color='grey', ls='', marker='o', mew=0, ms=5)
        
    ax.set_aspect('equal')
    if opts['GEOMETRY'] == 'cylindrical':
        ax.set_xlim(-rmax, rmax)
        ax.set_ylim(-rmax, rmax)
    else:
        ax.set_xlim(rjph.min(), rjph.max())
        ax.set_ylim(pimh_min, piph.max())
    
    if colorbar:
        cb = fig.colorbar(C)
        cb.set_label(label, fontsize=24)

    if xlabel == None:
        xlabel = r'$x$'
    if ylabel == None:
        ylabel = r'$y$'

    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)

def plotPhiSlice(fig, ax, rjph, zkph, q, label, pars, vmin=None, vmax=None, 
            noGhost=False, colorbar=True, xlabel=None, ylabel=None, log=False, 
            rmax=None, planets=None, cmap=None):

    if cmap is None:
        try:
            cmap = mpl.cm.inferno
        except:
            cmap = mpl.cm.afmhot
        
    if vmin is None:
        vmin = q.min()
    if vmax is None:
        vmax = q.max()

    Nr = rjph.shape[0]-1
    Nz = zkph.shape[0]-1

    if noGhost:
        if rjph[0] == 0:
            Nr_loc = Nr - 2
            j1 = 0
        else:
            Nr_loc = Nr - 4
            j1 = 2
        if Nz_loc < 1:
            Nz_loc = 1
            k1 = 0
        else:
            Nz_loc = Nz - 4
            k1 = 2
        dat = np.empty((Nz_loc,Nr_loc))
        for k in range(Nz_loc):
            i = Nr*(k+k1) + j1
            dat[k,:] = q[k+k1,j1:j1+Nr_loc]
        rjph_loc = rjph[j1:j1+Nr_loc+1]
        zkph_loc = zkph[k1:k1+Nz_loc+1]
    else:
        dat = q.copy()
        rjph_loc = rjph.copy()
        zkph_loc = zkph.copy()


    if rmax is None or rmax <=0.0:
        lim_float = True
    else:
        lim_float = False

    if log:
        norm = mpl.colors.LogNorm(vmin, vmax)
    else:
        norm = mpl.colors.Normalize(vmin, vmax)
        
    C = ax.pcolormesh(rjph_loc, zkph_loc, dat, 
            cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)

    if lim_float and rjph_loc.max() > rmax:
        rmax = rjph_loc.max()

    if planets is not None:
        rpl = planets[:,3]
        ppl = planets[:,4]
        xpl = rpl * np.cos(ppl)
        ypl = rpl * np.sin(ppl)
        ax.plot(xpl, ypl, color='grey', ls='', marker='o', mew=0, ms=5)
        
    ax.set_aspect('equal')
    ax.set_xlim(rjph.min(), rmax)
    ax.set_ylim(zkph.min(), zkph.max())
    
    if colorbar:
        cb = fig.colorbar(C)
        cb.set_label(label, fontsize=24)

    if xlabel == None:
        xlabel = r'$x$'
    if ylabel == None:
        ylabel = r'$y$'

    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)

def calcBounds(files):

    print("Checking bounds...")
    f = files[0]
    
    print("  Checking bounds for {0:s}...".format(f))
    prim = du.loadCheckpointPrims(f)

    num_q = prim.shape[1]
    bounds = [[prim[:,q].min(), prim[:,q].max()] for q in range(num_q)]
    bounds = np.array(bounds)

    for f in files[1:]:
        print("  Checking bounds for {0:s}...".format(f))
        prim = du.loadCheckpointPrims(f)
        for q in range(num_q):
            bounds[q,0] = min(bounds[q,0], prim[:,q].min())
            bounds[q,1] = max(bounds[q,1], prim[:,q].max())
    
    return bounds

def writeBoundsFile(filename, names, bounds):

    lines = ["{0:s} {1:.12g} {2:.12g}\n".format(names[q], bounds[q,0], bounds[q,1])
                for q in range(len(names))]

    f = open(filename, "w")
    for line in lines:
        f.write(line)

    f.close()

def readBoundsFile(filename, num_q):

    bounds = np.loadtxt(filename, usecols=[1,2])

    bounds = bounds[:num_q]

    return bounds
    
