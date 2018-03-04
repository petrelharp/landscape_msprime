import msprime
import numpy as np
import itertools
import random

def samples_without_replacement(x, k):
    combns = itertools.combinations(x, k)
    for c in combns:
        for d in itertools.permutations(c):
            yield d


def random_filtered_samples(f, x, k, n, replace=True, max_resamples=1e3):
    """
    Returns n randomly chosen *nonrepeating* samples of size k
    from x, sampling with or without replacement, filtered by f.
    If n is None, returns all of them.
    """
    if replace and (k > len(x)):
        raise ValueError("Not enough elements of x.")
    if replace:
        num_samples = len(x)**k
    else:
        num_samples = np.prod([len(x) - j for j in range(k)])
    if n is None:
        n = num_samples
    if n > num_samples:
        raise ValueError("There aren't that many samples.")
    if (num_samples < 1e6) or (n > 0.125 * num_samples):
        if replace:
            out = list(filter(f, itertools.product(x, repeat=k)))
        else:
            out = list(filter(f, samples_without_replacement(x, k)))
        random.shuffle(out)
        return out[:n]
    else:
        n_resamples = 0
        out = []
        j = 0
        while len(out) < n and n_resamples < max_resamples:
            samples = [np.random.choice(x, size=k, replace=replace) 
                       for _ in range(2*n)]
            out += list(filter(f, samples))
            n_resamples += 1
        # this should raise some actual error
        assert(n_resamples < max_resamples)
        return out[:n]


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
