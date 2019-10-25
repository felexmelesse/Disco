import sys
import h5py as h5
import discopy.util as util
import numpy as np


def diffOpts(filename1, filename2):

    opts1 = util.loadOpts(filename1)
    opts2 = util.loadOpts(filename2)

    if opts1 != opts2:
        print("Opts differences")


def diffPars(filename1, filename2):

    pars1 = util.loadPars(filename1)
    pars2 = util.loadPars(filename2)

    if pars1 != pars2:
        print("Pars differences")


def diffGit(filename1, filename2):

    hash1 = util.loadGitVersion(filename1)
    hash2 = util.loadGitVersion(filename2)

    if hash1 != hash2:
        print("Git Version differences")


def diffGrid(filename1, filename2):

    f1 = h5.File(filename1, "r")
    Np1 = f1['Grid/Np'][...]
    T1 = f1['Grid/T'][0]
    rjph1 = f1['Grid/r_jph'][...]
    zkph1 = f1['Grid/z_kph'][...]
    f1.close()

    f2 = h5.File(filename2, "r")
    Np2 = f2['Grid/Np'][...]
    T2 = f2['Grid/T'][0]
    rjph2 = f2['Grid/r_jph'][...]
    zkph2 = f2['Grid/z_kph'][...]
    f2.close()

    if T1 != T2:
        print("T difference")

    if Np1.shape != Np2.shape or (Np1 != Np2).any():
        print("Np differences")

    if rjph1.shape != rjph2.shape or (rjph1 != rjph2).any():
        print("r_jph differences")

    if zkph1.shape != zkph2.shape or (zkph1 != zkph2).any():
        print("z_kph differences")


def diffData(filename1, filename2):

    f1 = h5.File(filename1, "r")
    idPhi01 = f1['Grid/Id_phi0'][...]
    index1 = f1['Grid/Index'][...]
    Np1 = f1['Grid/Np'][...]
    cells1 = f1['Data/Cells'][...]
    diag1 = f1['Data/Diagnostics'][...]
    pl1 = f1['Data/Planets'][...]
    f1.close()

    f2 = h5.File(filename2, "r")
    idPhi02 = f2['Grid/Id_phi0'][...]
    index2 = f2['Grid/Index'][...]
    Np2 = f2['Grid/Np'][...]
    cells2 = f2['Data/Cells'][...]
    diag2 = f2['Data/Diagnostics'][...]
    pl2 = f2['Data/Planets'][...]
    f2.close()

    Nz = index1.shape[0]
    Nr = index1.shape[1]

    good = True

    for k in range(Nz):
        for j in range(Nr):
            a1 = index1[k, j]
            b1 = index1[k, j] + Np1[k, j]
            ann1 = cells1[a1:b1, :]
            a2 = index2[k, j]
            b2 = index2[k, j] + Np2[k, j]
            ann2 = cells2[a2:b2, :]
            if ann1.shape != ann2.shape or (ann1 != ann2).any():
                good = False
                break

    if not good:
        print("Cells differences")

    ph01 = np.empty((Nz, Nr, cells1.shape[1]))
    ph02 = np.empty((Nz, Nr, cells2.shape[1]))

    for k in range(Nz):
        for j in range(Nr):
            ph01[k, j, :] = cells1[idPhi01[k, j]]
            ph02[k, j, :] = cells2[idPhi02[k, j]]

    if ph01.shape != ph02.shape or (ph01 != ph02).any():
        print("Phi0 differences")

    if diag1.shape != diag2.shape or (diag1 != diag2).any():
        print("Diagnostics differences")

    if pl1.shape != pl2.shape or (pl1 != pl2).any():
        print("Plants differences")


def diff(filename1, filename2):

    diffGit(filename1, filename2)
    diffOpts(filename1, filename2)
    diffPars(filename1, filename2)
    diffGrid(filename1, filename2)
    diffData(filename1, filename2)


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Need 2 checkpoints bub!")
        sys.exit()

    filename1 = sys.argv[1]
    filename2 = sys.argv[2]

    diff(filename1, filename2)
