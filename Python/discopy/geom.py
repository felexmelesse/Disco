import numpy as np


def getXYZ(x1, x2, x3, opts, pars):

    x1, x2, x3 = np.broadcast_arrays(x1, x2, x3)

    if opts['GEOMETRY'] == "cylindrical":
        return getXYZcyl(x1, x2, x3, pars)
    elif opts['GEOMETRY'] == "spherical":
        return getXYZsph(x1, x2, x3, pars)
    else:
        return getXYZcart(x1, x2, x3, pars)


def getCellCentroids(dat, opts, pars):
    x1jph = dat[0]
    x2iph = dat[3]
    x3kph = dat[1]
    nx2 = dat[7]
    index = dat[8]
    x2max = pars['Phi_Max']
    N = x2iph.shape[0]

    nx1 = nx2.shape[1]
    nx3 = nx2.shape[0]

    x11d = getCentroid(x1jph[:-1], x1jph[1:], 1, opts)
    x31d = getCentroid(x3kph[:-1], x3kph[1:], 3, opts)

    x1 = np.empty(N)
    x2 = np.empty(N)
    x3 = np.empty(N)

    for k in range(nx3):
        for j in range(nx1):
            a = index[k, j]
            b = index[k, j] + nx2[k, j]

            x1[a:b] = x11d[j]
            x3[a:b] = x31d[k]

            x2p = x2iph[a:b].copy()
            x2p[x2p < 0.0] += x2max
            x2p[x2p >= x2max] -= x2max
            x2m = np.roll(x2p, 1)
            x2m[x2m > x2p] -= x2max

            x2[a:b] = getCentroid(x2m, x2p, 2, opts)

    return x1, x2, x3


def getCentroid(xm, xp, dim, opts):
    if opts['GEOMETRY'] == "cylindrical":
        return getCentroidcyl(xm, xp, dim)
    elif opts['GEOMETRY'] == "spherical":
        return getCentroidsph(xm, xp, dim)
    else:
        return getCentroidcart(xm, xp, dim)


def getDV(dat, opts, pars):

    x1jph = dat[0]
    x2iph = dat[3]
    x3kph = dat[1]
    nx2 = dat[7]
    index = dat[8]

    if opts['GEOMETRY'] == "cylindrical":
        return getDVcyl(x1jph, x2iph, x3kph, nx2, index, pars)
    elif opts['GEOMETRY'] == "spherical":
        return getDVsph(x1jph, x2iph, x3kph, nx2, index, pars)
    else:
        return getDVcart(x1jph, x2iph, x3kph, nx2, index, pars)


def getVXYZ(x1, x2, x3, v1, v2, v3, opts):

    h1, h2, h3 = getScaleFactors(x1, x2, x3, opts)
    return getVecXYZ(x1, x2, x3, h1*v1, h2*v2, h3*v3, opts)


def getBXYZ(x1, x2, x3, B1, B2, B3, opts):

    return getVecXYZ(x1, x2, x3, B1, B2, B3, opts)


def getVecXYZ(x1, x2, x3, v1, v2, v3, opts):

    if opts['GEOMETRY'] == "cylindrical":
        return getVecXYZcyl(x1, x2, x3, v1, v2, v3)
    else:
        return getVecXYZcart(x1, x2, x3, v1, v2, v3)


def getScaleFactors(x1, x2, x3, opts):
    if opts['GEOMETRY'] == "cylindrical":
        return getScaleFactorscyl(x1, x2, x3)
    else:
        return getScaleFactorscart(x1, x2, x3)


def getXYZcart(x, y, z, pars):
    return x, y, z


def getXYZcyl(r, phi, z, pars):
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    return x, y, z


def getXYZsph(r, phi, th, pars):
    x = r*np.cos(phi)*np.sin(th)
    y = r*np.sin(phi)*np.sin(th)
    z = r*np.cos(th)
    return x, y, z


def getCentroidcart(xm, xp, dim):
    return 0.5*(xm+xp)


def getCentroidcyl(xm, xp, dim):
    if dim == 1:
        xc = 2.0*(xm*xm + xm*xp + xp*xp) / (3.0*(xm+xp))
    else:
        xc = 0.5*(xm+xp)
    return xc


def getCentroidsph(xm, xp, dim):
    if dim == 1:
        xc = 3.0*(xm*xm*xm + xm*xm*xp + xm*xp*xp + xp*xp*xp
                  ) / (4.0*(xm*xm+xm*xp+xp*xp))
    elif dim == 3:
        xc = np.empty(xm.shape)
        xc[xm == xp] = xm[xm == xp]
        dif = xm != xp
        cp = np.cos(xp)
        cm = np.cos(xm)
        sp = np.sin(xp)
        sm = np.sin(xm)
        xc[dif] = ((xm*cm-xp*cp+sp-sm) / (cm-cp))[dif]
    else:
        xc = 0.5*(xm+xp)
    return xc


def getDVcart(xjph, yiph, zkph, ny, index, pars):

    nz = ny.shape[0]
    nx = ny.shape[1]

    dx1d = xjph[1:] - xjph[:-1]
    dz1d = zkph[1:] - zkph[:-1]

    ymax = pars['Phi_Max']

    dx = np.empty(yiph.shape)
    dy = np.empty(yiph.shape)
    dz = np.empty(yiph.shape)

    for k in range(nz):
        for j in range(nx):
            a = index[k, j]
            b = index[k, j]+ny[k, j]

            dx[a:b] = dx1d[j]
            dz[a:b] = dz1d[k]

            yp = yiph[a:b]
            ym = np.roll(yp, 1)
            ym[ym > yp] -= ymax
            dy[a:b] = yp-ym

    if nz > 1:
        dV = dx*dy*dz
    else:
        dV = dx*dy

    return dV


def getDVcyl(rjph, piph, zkph, nphi, index, pars):

    nz = nphi.shape[0]
    nr = nphi.shape[1]

    r1d = 0.5*(rjph[1:] + rjph[:-1])
    dr1d = rjph[1:] - rjph[:-1]
    dz1d = zkph[1:] - zkph[:-1]

    pmax = pars['Phi_Max']

    r = np.empty(piph.shape)
    dr = np.empty(piph.shape)
    dp = np.empty(piph.shape)
    dz = np.empty(piph.shape)

    for k in range(nz):
        for j in range(nr):
            a = index[k, j]
            b = index[k, j]+nphi[k, j]

            r[a:b] = r1d[j]
            dr[a:b] = dr1d[j]
            dz[a:b] = dz1d[k]

            pp = piph[a:b]
            pm = np.roll(pp, 1)
            pm[pm > pp] -= pmax
            dp[a:b] = pp-pm

    if nz > 1:
        dV = r*dr*dp*dz
    else:
        dV = r*dr*dp

    return dV


def getDVsph(rjph, piph, tkph, nphi, index, pars):

    nt = nphi.shape[0]
    nr = nphi.shape[1]

    r21d = (rjph[1:]*rjph[1:] + rjph[1:]*rjph[:-1] + rjph[:-1]*rjph[:-1])/3.0
    sinth1d = np.sin(0.5*(tkph[:-1]+tkph[1:]))
    dr1d = rjph[1:] - rjph[:-1]
    sindth1d = 2.0*np.sin(0.5*(tkph[1:] - tkph[:-1]))

    pmax = pars['Phi_Max']

    r2 = np.empty(piph.shape)
    sinth = np.empty(piph.shape)
    dr = np.empty(piph.shape)
    dp = np.empty(piph.shape)
    sindth = np.empty(piph.shape)

    for k in range(nt):
        for j in range(nr):
            a = index[k, j]
            b = index[k, j]+nphi[k, j]

            r2[a:b] = r21d[j]
            dr[a:b] = dr1d[j]
            sinth[a:b] = sinth1d[k]
            sindth[a:b] = sindth1d[k]

            pp = piph[a:b]
            pm = np.roll(pp, 1)
            pm[pm > pp] -= pmax
            dp[a:b] = pp-pm

    dV = r2*sinth*dr*dp*sindth

    return dV


def integrate(f, dat, opts, pars, dV=None):
    if dV is None:
        dV = getDV(dat, opts, pars)

    nphi = dat[7]
    index = dat[8]

    mask = np.empty(f.shape, dtype=np.bool)
    mask[:] = True

    nz = nphi.shape[0]
    nr = nphi.shape[1]

    for k in range(nz):
        for j in range(nr):
            a = index[k, j]
            b = index[k, j] + nphi[k, j]
            if (j < 2 and pars["NoBC_Rmin"] == 0)\
                    or (j >= nr-2 and pars["NoBC_Rmax"] == 0)\
                    or (k < 2 and nz > 1 and pars["NoBC_Zmin"] == 0)\
                    or (k >= nz-2 and nz > 1 and pars["NoBC_Zmax"] == 0):
                mask[a:b] = False

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


def getSignedDphi(phip, phim, pars):

    pmax = pars['Phi_Max']

    dphi = phip-phim
    while dphi > 0.5*pmax:
        dphi -= pmax
    while dphi < -0.5*pmax:
        dphi += pmax

    return dphi


def calcDivB(dat, opts, pars, dV=None):
    if dV is None:
        dV = getDV(dat, opts, pars)

    Phi = dat[5]
    nphi = dat[7]
    index = dat[8]
    piph = dat[3]
    idPhi0 = dat[6]

    nz = nphi.shape[0]
    nr = nphi.shape[1]

    divB = np.zeros(piph.shape)

    PhiP = Phi[:, 0]
    PhiRm = Phi[:, 1]
    PhiRp = Phi[:, 2]
    PhiZm = Phi[:, 3]
    PhiZp = Phi[:, 4]

    # phi direction
    for k in range(nz):
        for j in range(nr):
            a = index[k, j]
            b = index[k, j] + nphi[k, j]

            divB[a:b] += PhiP[a:b] - np.roll(PhiP[a:b], 1)

    # r direction
    for k in range(nz):
        for jL in range(nr-1):
            jR = jL+1
            aL = index[k, jL]
            bL = index[k, jL] + nphi[k, jL]
            aR = index[k, jR]
            bR = index[k, jR] + nphi[k, jR]
            piphL = piph[aL:bL]
            piphR = piph[aR:bR]
            PhiRL = PhiRp[aL:bL]
            PhiRR = PhiRm[aR:bR]

            iL0 = idPhi0[k, jL] - aL
            iR0 = idPhi0[k, jR] - aR
            iL = iL0
            iR = iR0
            done = False

            while not done:
                phiL = piphL[iL]
                phiR = piphR[iR]

                dp = getSignedDphi(phiL, phiR, pars)

                if dp < 0.0:
                    phi_back = np.roll(piphR, 1)
                    dp1 = getSignedDphi(phiL, phi_back[iR], pars)
                    if dp1 <= 0.0:
                        print("Uh Oh")
                if dp > 0.0:
                    phi_back = np.roll(piphL, 1)
                    dp1 = getSignedDphi(phiR, phi_back[iL], pars)
                    if dp1 <= 0.0:
                        print("Oh Uh")

                if dp == 0.0:
                    divB[aL+iL] += PhiRL[iL]
                    divB[aR+iR] -= PhiRL[iL]
                    iL += 1
                    iR += 1
                elif dp < 0.0:
                    divB[aL+iL] += PhiRL[iL]
                    divB[aR+iR] -= PhiRL[iL]
                    iL += 1
                else:
                    divB[aL+iL] += PhiRR[iR]
                    divB[aR+iR] -= PhiRR[iR]
                    iR += 1

                if iL == nphi[k, jL]:
                    iL = 0
                if iR == nphi[k, jR]:
                    iR = 0

                if iL == iL0 and iR == iR0:
                    done = True
    # z direction

    for kL in range(nz-1):
        for j in range(nr):
            kR = kL+1
            aL = index[kL, j]
            bL = index[kL, j] + nphi[kL, j]
            aR = index[kR, j]
            bR = index[kR, j] + nphi[kR, j]
            piphL = piph[aL:bL]
            piphR = piph[aR:bR]
            PhiZL = PhiZp[aL:bL]
            PhiZR = PhiZm[aR:bR]

            pimhL = np.roll(piphL, 1)
            pimhR = np.roll(piphR, 1)

            iL0 = idPhi0[kL, j] - aL
            iR0 = idPhi0[kR, j] - aR
            iL = iL0
            iR = iR0
            done = False

            while not done:
                phiL = piphL[iL]
                phiR = piphR[iR]

                dp = getSignedDphi(phiL, phiR, pars)

                if dp < 0.0:
                    dp1 = getSignedDphi(phiL, pimhR[iR], pars)
                    if dp1 <= 0.0:
                        print("Uh Oh")
                if dp > 0.0:
                    dp1 = getSignedDphi(phiR, pimhL[iL], pars)
                    if dp1 <= 0.0:
                        print("Oh Uh")

                if dp == 0.0:
                    print("quad?")
                    divB[aL+iL] += PhiZL[iL]
                    divB[aR+iR] -= PhiZL[iL]
                    iL += 1
                    iR += 1
                elif dp < 0.0:
                    divB[aL+iL] += PhiZL[iL]
                    divB[aR+iR] -= PhiZL[iL]
                    iL += 1
                else:
                    divB[aL+iL] += PhiZR[iR]
                    divB[aR+iR] -= PhiZR[iR]
                    iR += 1

                if iL == nphi[kL, j]:
                    iL = 0
                if iR == nphi[kR, j]:
                    iR = 0

                if iL == iL0 and iR == iR0:
                    done = True

    divB /= dV

    return divB


def getDA1(x1, x2ph, x3f, dat, opts, pars):

    x2max = pars['Phi_Max']
    index = dat[8]
    nx2 = dat[7]
    x2mh = np.empty(x2ph.shape)
    x3mh = np.empty(x2ph.shape)
    x3ph = np.empty(x2ph.shape)
    for k in range(index.shape[0]):
        for j in range(index.shape[1]):
            a = index[k, j]
            b = index[k, j] + nx2[k, j]
            x2mh[a:b] = np.roll(x2ph[a:b], 1)
            x3mh[a:b] = x3f[k]
            x3ph[a:b] = x3f[k+1]
    x2mh[x2mh > x2ph] -= x2max

    if opts['GEOMETRY'] == "cylindrical":
        return getDA1cyl(x1, x2ph, x2mh, x3ph, x3mh)
    elif opts['GEOMETRY'] == "spherical":
        return getDA1sph(x1, x2ph, x2mh, x3ph, x3mh)
    else:
        return getDA1cart(x1, x2ph, x2mh, x3ph, x3mh)


def getDA1cart(x1, x2ph, x2mh, x3ph, x3mh):
    return (x2ph-x2mh)*(x3ph-x3mh)


def getDA1cyl(x1, x2ph, x2mh, x3ph, x3mh):
    return x1*(x2ph-x2mh)*(x3ph-x3mh)


def getDA1sph(x1, x2ph, x2mh, x3ph, x3mh):
    sinth = np.sin(0.5*(x3ph+x3mh))
    sindth = 2*np.sin(0.5*(x3ph-x3mh))
    return x1*x1*sinth*sindth*(x2ph-x2mh)


def integrate2(f, dat, opts, pars):

    index = dat[8]
    n1 = index.shape[1]
    n2 = dat[7]
    n3 = index.shape[0]
    x2ph = dat[3]
    x2max = pars['Phi_Max']

    integral = np.empty((n3, n1), dtype=f.dtype)

    for k in range(n3):
        for j in range(n1):
            a = index[k, j]
            b = index[k, j] + n2[k, j]
            x2p = x2ph[a:b]
            x2m = np.roll(x2p, 1)
            x2m[x2m > x2p] -= x2max
            dx2 = x2p-x2m
            integral[k, j] = (f[a:b]*dx2).sum()

    if n3 == 1:
        integral = integral[0, :]

    return integral


def integrateTrans1(x1, f, dat, opts, pars, dA=None):

    if dA is None:
        dA = getDA1(x1, dat[3], dat[1], dat, opts, pars)

    mask = np.empty(f.shape, dtype=np.bool)
    mask[:] = True
    x1ph = dat[0]
    index = dat[8]
    n1 = index.shape[1]
    n2 = dat[7]
    n3 = index.shape[0]

    for k in range(n3):
        for j in range(n1):
            a = index[k, j]
            b = index[k, j] + n2[k, j]
            if (k < 2 and n3 > 1 and pars["NoBC_Zmin"] == 0)\
                    or (k >= n3-2 and n3 > 1 and pars["NoBC_Zmax"] == 0):
                mask[a:b] = False

    integral = np.empty(n1, dtype=f.dtype)
    for j in range(n1):
        ind = (x1 < x1ph[j+1]) & (x1 > x1ph[j])
        integral[j] = (f[mask*ind]*dA[mask*ind]).sum()

    return integral