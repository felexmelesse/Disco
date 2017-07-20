import sys
import math
import h5py as h5
import numpy as np

strType = h5.special_dtype(vlen=bytes)

def loadPars(f):
    pars = dict()

    parG = f['Pars']

    for key in parG:
        pars[key] = parG[key][0]

    return pars

def loadOpts(f):
    opts = dict()

    optG = f['Opts']

    for key in optG:
        opts[key] = optG[key][0]

    return opts

def loadGrid(f):
    g = f['Grid']
    id_phi0 = g['Id_phi0'][...]
    index = g['Index'][...]
    nphi = g['Np'][...]
    t = g['T'][0]
    r_jph = g['r_jph'][...]
    z_kph = g['z_kph'][...]

    grid = (t, id_phi0, index, nphi, r_jph, z_kph)

    return grid

def loadData(f):

    g = f['Data']
    ct_mode = f['Opts']['CT_MODE'][0]
    if ct_mode == 1:
        nfaces = 3
        bflux = g['Cells'][:,-1-nfaces:-1]
    elif ct_mode == 2:
        nfaces = 5
        bflux = g['Cells'][:,-1-nfaces:-1]
    else:
        nfaces = 0
        bflux = None

    prim = g['Cells'][:,:-1-nfaces]
    piph = g['Cells'][:,-1]
    planets = g['Planets'][...]
    
    try:
        diagRZ = g['Poloidal_Diagnostics'][...]
    except KeyError:
        diagRZ = None
    try:
        diagR = g['Radial_Diagnostics'][...]
    except KeyError:
        diagR = None
    try:
        diagZ = g['Vertical_Diagnostics'][...]
    except KeyError:
        diagZ = None

    data = (prim, bflux, piph, planets, diagRZ, diagR, diagZ)

    return data

def savePars(f, pars):

    g = f.create_group("Pars")
    for key in pars:
        g.create_dataset(key, data=np.array([pars[key]]))

def saveOpts(f, opts):

    g = f.create_group("Opts")
    for key in opts:
        if isinstance(opts[key], str):
            d = g.create_dataset(key, (1,), dtype=strType)
            d[0] = opts[key]
        else:
            g.create_dataset(key, data=np.array([opts[key]]))

def saveGrid(f, grid):

    g = f.create_group("Grid")
    g.create_dataset('T', data=np.array([grid[0]]))
    g.create_dataset('Id_phi0', data=grid[1])
    g.create_dataset('Index', data=grid[2])
    g.create_dataset('Np', data=grid[3])
    g.create_dataset('r_jph', data=grid[4])
    g.create_dataset('z_kph', data=grid[5])

def saveData(f, data):

    g = f.create_group("Data")
    g.create_dataset('Planets', data=data[3])
    if data[4] is not None:
        g.create_dataset('Poloidal_Diagnostics', data=data[4])
    if data[5] is not None:
        g.create_dataset('Radial_Diagnostics', data=data[5])
    if data[6] is not None:
        g.create_dataset('Vertical_Diagnostics', data=data[6])

    prim = data[0]
    bflux = data[1]
    piph = data[2]

    if bflux is None:
        nfaces = 0
    else:
        nfaces = bflux.shape[1]
    nq = prim.shape[1]
    N = prim.shape[0]

    d = g.create_dataset('Cells', (N,nq+nfaces+1), dtype=np.float64)
    d[:,:nq] = prim[:,:]
    if bflux is not None:
        d[:,nq:nq+nfaces] = bflux[:,:]
    d[:,-1] = piph[:]


def loadCheckpointAll(filename):

    f = h5.File(filename, "r")

    GIT_VERSION = f['GIT_VERSION'][0]
    pars = loadPars(f)
    opts = loadOpts(f)
    grid = loadGrid(f)
    data = loadData(f)

    f.close()

    checkpoint = (GIT_VERSION, pars, opts, grid, data)

    return checkpoint

def loadCheckpointData(filename):

    f = h5.File(filename, "r")

    piph = f['Data']['Cells'][:,-1][...]
    prim = f['Data']['Cells'][:,:-1][...]
    index = f['Grid']['Index'][...]
    idPhi0 = f['Grid']['Id_phi0'][...]
    nphi = f['Grid']['Np'][...]
    t = f['Grid']['T'][0]
    riph = f['Grid']['r_jph'][...]
    ziph = f['Grid']['z_kph'][...]

    r = np.zeros(piph.shape)
    z = np.zeros(piph.shape)
    phi = np.zeros(piph.shape)
    R = 0.5*(riph[1:] + riph[:-1])
    Z = 0.5*(ziph[1:] + ziph[:-1])
    primPhi0 = np.zeros((index.shape[0], index.shape[1], prim.shape[1]))
    for k in range(index.shape[0]):
        for j in range(index.shape[1]):
            ind0 = index[k,j]
            ind1 = ind0 + nphi[k,j]
            r[ind0:ind1] = R[j]
            z[ind0:ind1] = Z[k]
            piph_strip = piph[ind0:ind1]
            pimh = np.roll(piph_strip, 1)
            pimh[pimh>piph_strip] -= 2*np.pi
            phi[ind0:ind1] = 0.5*(pimh+piph_strip)
            primPhi0[k,j,:] = prim[idPhi0[k,j],:]

    return t, r, phi, z, prim, (riph, ziph, primPhi0, piph)

def loadDiagRZ(filename):

    f = h5.File(filename, "r")

    t = f['Grid']['T'][0]
    rjph = f['Grid']['r_jph'][...]
    zkph = f['Grid']['z_kph'][...]
    diag = f['Data']['Poloidal_Diagnostics'][...]

    f.close()

    Nr = rjph.shape[0]-1
    Nz = zkph.shape[0]-1
    Nq = diag.shape[-1]
    #diag = np.resize(diag, (Nz,Nr,Nq))

    R = 0.5*(rjph[1:]+rjph[:-1])
    Z = 0.5*(zkph[1:]+zkph[:-1])

    r = np.empty((Nz,Nr))
    z = np.empty((Nz,Nr))

    r[:,:] = R[None,:]
    z[:,:] = Z[:,None]

    return t, r, z, diag, rjph, zkph

def saveCheckpoint(checkpoint, filename):

    git = checkpoint[0]
    pars = checkpoint[1]
    opts = checkpoint[2]
    grid = checkpoint[3]
    data = checkpoint[4]

    f = h5.File(filename, "w")

    d = f.create_dataset("GIT_VERSION", (1,), dtype=strType)
    d[0] = git

    savePars(f, pars)
    saveOpts(f, opts)
    saveGrid(f, grid)
    saveData(f, data)

    f.close()

def setTime(filename, t=0.0):
    f = h5.File(filename, "r+")
    f['Grid/T'][0] = t
    f.close()

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: $ python checkpoint.py <checkpoint.h5>")
        print("Loads checkpoint data and prints to screen.")
        sys.exit()

    checkpoint = loadCheckpointAll(sys.argv[1])

    print("GIT")
    print(checkpoint[0])
    print("PARS")
    print(checkpoint[1])
    print("OPTS")
    print(checkpoint[2])
    print("GRID")
    print(checkpoint[3])
    print("DATA")
    print(checkpoint[4])

    saveCheckpoint(checkpoint, "discopytest.h5")

    checkpoint2 = loadCheckpointAll("discopytest.h5")

    print("GIT")
    print(checkpoint2[0])
    print("PARS")
    print(checkpoint2[1])
    print("OPTS")
    print(checkpoint2[2])
    print("GRID")
    print(checkpoint2[3])
    print("DATA")
    print(checkpoint2[4])

