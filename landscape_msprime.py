import msprime
import numpy as np
import itertools
import random

from tools import *

def read_migration_matrix(fname):
    f = open(fname, "r")
    header = f.readline().split()
    f.close()
    i = np.loadtxt(fname, skiprows=1, usecols=header.index('"i"'), dtype='int16')
    j = np.loadtxt(fname, skiprows=1, usecols=header.index('"j"'), dtype='int16')
    x = np.loadtxt(fname, skiprows=1, usecols=header.index('"x"'))
    ids = list(set(list(set(i)) + list(set(j))))
    ids.sort()
    M = np.empty((len(ids), len(ids)))
    for ii, jj, xx in zip(i, j, x):
        M[ii, jj] = xx
    # msprime wants diag elements to be zero
    for ii in ids:
        M[ii, ii] = 0.0
    return M

def generate_stat_indices(stat, num_sets, nstats=None):
    """
    Given the name of a statistic and a number of sample_sets,
    generate a list of indices suitable to be passed into a
    statistics calculator. If nstats is None, then will find *all*
    nonredundant combinations; otherwise will generate a random
    collection of these.

    Possible statistics:

    - f2: (a,b) = (b,a)
    - f3: (a;b,c) = (a;c,b)
    - f4: (a,b;c,d) = (c,d;a,b) = -(b,a;c,d)
    - Y2: (a,b)
    - Y3: (a;b,c) = (a;c,b)

    """
    if stat == "f2":
        indices = random_filtered_samples(lambda x: x[0] < x[1], 
                                          x=range(num_sets), k=2,
                                          n=nstats, replace=False)
    elif stat == "f3":
        indices = random_filtered_samples(lambda x: x[1] < x[2], 
                                          x=range(num_sets), k=3,
                                          n=nstats, replace=False)
    elif stat == "f4":
        indices = random_filtered_samples(lambda x: (x[0] < x[1]) 
                                                     & (x[2] < x[3]) 
                                                     & (x[0] < x[2]), 
                                          x=range(num_sets), k=4,
                                          n=nstats, replace=False)
    elif stat == "Y2":
        indices = random_filtered_samples(lambda x: True,
                                          x=range(num_sets), k=2,
                                          n=nstats, replace=False)
    elif stat == "Y3":
        indices = random_filtered_samples(lambda x: x[1] < x[2], 
                                          x=range(num_sets), k=3,
                                          n=nstats, replace=False)
    else:
        raise ValueError("Unknown statistic" + stat)
    return indices


def compute_stat(ts, statname, nstats):
    pops = { u : ts.node(u).population for u in ts.samples()}
    unique_pops = list(set(pops.values()))
    unique_pops.sort()
    A = [[u for u in pops if pops[u] == a] for a in unique_pops]
    branch_stats = msprime.stats.BranchLengthStatCalculator(ts)
    site_stats = msprime.stats.SiteStatCalculator(ts)
    branch_fn = getattr(branch_stats, statname + "_vector")
    site_fn = getattr(site_stats, statname + "_vector")
    print("generating indices")
    inds = generate_stat_indices(statname, num_sets=len(unique_pops), 
                                 nstats=nstats)
    x = np.empty((nstats, 2))
    print("computing stats")
    x[:,0] = branch_fn(A, [0, ts.sequence_length], inds)[0]
    x[:,1] = site_fn(A, [0, ts.sequence_length], inds)[0]
    return x


def f(x) : 
    return (x[0] < x[1]) & (x[2] < x[3]) & (x[0] < x[2])

