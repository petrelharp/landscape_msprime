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



