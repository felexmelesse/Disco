import math
import numpy as np


def getXYZ(x1, x2, x3, opts, pars):

    x1, x2, x3 = np.broadcast_arrays(x1, x2, x3)

    if opts['GEOMETRY'] == "Cylindrical":
        return getXYZcyl(x1, x2, x3, pars)
    else:
        return getXYZcart(x1, x2, x3, pars)

def getCellCentroids(dat, opts, pars):
    x1jph = dat[0]
    x2iph = dat[3]
    x3kph = dat[1]
    nx2 = dat[7]
    x2max = pars['Phi_Max']
    N = x2iph.shape[0]

    nx1 = nx2.shape[1]
    nx3 = nx2.shape[0]

    x11d = getCentroid(x1jph[:-1], x1jph[1:], 1, opts)
    x31d = getCentroid(x3kph[:-1], x3kph[1:], 3, opts)

    x1 = np.empty(N)
    x2 = np.empty(N)
    x3 = np.empty(N)

    i = 0
    for k in range(nx3):
        for j in range(nx1):
            a = i
            b = i + nx2[k,j]

            x1[a:b] = x11d[j]
            x3[a:b] = x31d[k]

            x2p = x2iph[a:b].copy()
            x2p[x2p < 0.0] += x2max
            x2p[x2p >= x2max] -= x2max
            x2m = np.roll(x2p, 1)
            x2m[x2m > x2p] -= x2max

            x2[a:b] = getCentroid(x2m, x2p, 2, opts)

            i += nx2[k,j]

    return x1, x2, x3

def getCentroid(xm, xp, dim, opts):
    if opts['GEOMETRY'] == "Cylindrical":
        return getCentroidcyl(xm, xp, dim)
    else:
        return getCentroidcart(xm, xp, dim)

def getDV(dat, opts, pars):

    x1jph = dat[0]
    x2iph = dat[3]
    x3kph = dat[1]
    nx2 = dat[7]

    if opts['GEOMETRY'] == "Cylindrical":
        return getDVcyl(x1jph, x2iph, x3kph, nx2, pars)
    else:
        return getDVcart(x1jph, x2iph, x3kph, nx2, pars)

def getVXYZ(x1, x2, x3, v1, v2, v3, opts):

    h1, h2, h3 = getScaleFactors(x1, x2, x3, opts)
    return getVecXYZ(x1, x2, x3, h1*v1, h2*v2, h3*v3, opts)

def getBXYZ(x1, x2, x3, B1, B2, B3, opts):

    return getVecXYZ(x1, x2, x3, B1, B2, B3, opts)

def getVecXYZ(x1, x2, x3, v1, v2, v3, opts):

    if opts['GEOMETRY'] == "Cylindrical":
        return getVecXYZcyl(x1, x2, x3, v1, v2, v3)
    else:
        return getVecXYZcart(x1, x2, x3, v1, v2, v3)


def getScaleFactors(x1, x2, x3, opts):
    if opts['GEOMETRY'] == "Cylindrical":
        return getScaleFactorscyl(x1, x2, x3)
    else:
        return getScaleFactorscart(x1, x2, x3)

def getXYZcart(x, y, z, pars):
    return x, y, z

def getXYZcyl(r, phi, z, pars):
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    return x, y, z

def getCentroidcart(xm, xp, dim):
    return 0.5*(xm+xp)

def getCentroidcyl(xm, xp, dim):
    if dim == 1:
        xc = 2.0*(xm*xm + xm*xp + xp*xp) / (3.0*(xm+xp))
    else:
        xc = 0.5*(xm+xp)
    return xc


def getDVcart(xjph, yiph, zkph, ny, pars):

    nz = ny.shape[0]
    nx = ny.shape[1]

    dx1d = xjph[1:] - xjph[:-1]
    dz1d = zkph[1:] - zkph[:-1]

    ymax = pars['Phi_Max']

    dx = np.empty(yiph.shape)
    dy = np.empty(yiph.shape)
    dz = np.empty(yiph.shape)

    i = 0
    for k in range(nz):
        for j in range(nx):
            a = i
            b = i+ny[k,j]

            dx[a:b] = dx1d[j]
            dz[a:b] = dz1d[k]

            yp = yiph[a:b]
            ym = np.roll(yp, 1)
            ym[ym>yp] -= ymax
            dy[a:b] = yp-ym

            i += ny[k,j]

    if nz > 1:
        dV = dx*dy*dz
    else:
        dV = dx*dy

    return dV

def getDVcyl(rjph, piph, zkph, np, pars):

    nz = np.shape[0]
    nr = np.shape[1]

    r1d = 0.5*(rjph[1:] + rjph[:-1])
    dr1d = rjph[1:] - rjph[:-1]
    dz1d = zkph[1:] - zkph[:-1]

    pmax = pars['Phi_Max']

    r = np.empty(piph.shape)
    dr = np.empty(piph.shape)
    dp = np.empty(piph.shape)
    dz = np.empty(piph.shape)

    i = 0
    for k in range(nz):
        for j in range(nr):
            a = i
            b = i+np[k,j]

            r[a:b] = r1d[j]
            dr[a:b] = dr1d[j]
            dz[a:b] = dz1d[k]

            pp = piph[a:b]
            pm = np.roll(pp, 1)
            pm[pm>pp] -= pmax
            dp[a:b] = pp-pm

            i += np[k,j]

    if nz > 1:
        dV = r*dr*dp*dz
    else:
        dV = r*dr*dp

    return dV

def integrate(f, dat, opts, pars, dV=None):
    if dV is None:
        dV = getDV(dat, opts, pars)

    nphi = dat[7]

    mask = np.empty(f.shape, dtype=np.bool)
    mask[:] = True

    nz = nphi.shape[0]
    nr = nphi.shape[1]



    i = 0
    for k in range(nz):
        for j in range(nr):
            a = i
            b = i + nphi[k,j]
            if (j < 2 and pars["NoBC_Rmin"] == 0)\
                    or (j >= nr-2 and pars["NoBC_Rmax"] == 0)\
                    or (k < 2 and nz > 1 and pars["NoBC_Zmin"] == 0)\
                    or (k >= nz-2 and nz > 1 and pars["NoBC_Zmax"] == 0):
                mask[a:b] = False
            i += nphi[k,j]

    return (f[mask]*dV[mask]).sum()

def getScaleFactorscart(x, y, z):

    hx = np.ones(y.shape)
    hy = np.ones(y.shape)
    hz = np.ones(y.shape)

    return hx, hy, hz

def getScaleFactorscyl(r, phi, z):

    hr = np.ones(phi.shape)
    hp = np.empty(phi.shape)
    hz = np.ones(phi.shape)

    hp[:] = r

    return hr, hp, hz

def getVecXYZcart(x, y, z, vx, vy, vz):
    return vx, vy, vz

def getVecXYZcyl(r, phi, z, vr, vp, vz):

    cp = np.cos(phi)
    sp = np.sin(phi)

    vx = cp*vr - sp*vp
    vy = sp*vr + cp*vp

    return vx, vy, vz
