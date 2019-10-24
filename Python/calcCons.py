import sys
import numpy as np
import discopy.util as util
import discopy.geom as geom

def calc(filename):

    t, r, phi, z,prim, dat = util.loadCheckpoint(filename)
    opts = util.loadOpts(filename)
    pars = util.loadPars(filename)

    gam = pars['Adiabatic_Index']

    h1, h2, h3 = geom.getScaleFactors(r, phi, z, opts)

    rho = prim[:,0]
    P = prim[:,1]
    u1 = prim[:,2]*h1
    u2 = prim[:,3]*h2
    u3 = prim[:,4]*h3

    s1 = rho*u1*h1
    s2 = rho*u2*h2
    s3 = rho*u3*h3

    v2 = u1*u1 + u2*u2 + u3*u3
    e = 0.5*rho*v2 + P/(gam-1)

    if opts['NUM_C'] >= 8:
        b1 = prim[:,5]
        b2 = prim[:,6]
        b3 = prim[:,7]
        B2 = b1*b1 + b2*b2 + b3*b3
        e += 0.5*B2

    M = geom.integrate(rho, dat, opts, pars)
    S1 = geom.integrate(s1, dat, opts, pars)
    S2 = geom.integrate(s2, dat, opts, pars)
    S3 = geom.integrate(s3, dat, opts, pars)
    E = geom.integrate(e, dat, opts, pars)

    print(M, S1, S2, S3, E)

    return np.array([M, S1, S2, S3, E])

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need some checkpoints!")
        sys.exit()

    val0 = calc(sys.argv[1])
    
    for f in sys.argv[2:-1]:
        calc(f)

    val1 = calc(sys.argv[-1])

    err = np.empty(val0.shape)

    nz = val0 != 0.0
    err[nz] = (val1[nz]-val0[nz]) / val0[nz] 
    err[~nz] = val1[~nz]

    print(err)
