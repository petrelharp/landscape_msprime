import numpy as np

def read_migration_matrix(fname):
    f = open(fname, "r")
    header = f.readline().split()
    f.close()
    i = np.loadtxt(fname, skiprows=1, usecols=header.index('"i"'), dtype='int16')
    j = np.loadtxt(fname, skiprows=1, usecols=header.index('"j"'), dtype='int16')
    x = np.loadtxt(fname, skiprows=1, usecols=header.index('"x"'))
    if (x < 0.0).any():
        raise ValueError("Migration rates must be positive.")
    ids = list(set(list(set(i)) + list(set(j))))
    ids.sort()
    M = np.empty((len(ids), len(ids)))
    for ii, jj, xx in zip(i, j, x):
        M[ii, jj] = xx
    # msprime wants diag elements to be zero
    for ii in ids:
        M[ii, ii] = 0.0
    return M
