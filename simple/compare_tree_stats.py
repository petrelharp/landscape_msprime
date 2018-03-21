import msprime
import numpy as np
import time

import landscape_msprime as lmsp


def run_sim(M, nsamples, **kwargs):
    assert(M.shape[0] == M.shape[1])
    n = M.shape[0]
    print("Simulating on a landscape of {} patches.".format(n))
    # sample a total of nsamples, uniformly spread
    sample_sizes = [int(np.ceil(nsamples/n)) for _ in range(n)]
    while sum(sample_sizes) < nsamples:
        sample_sizes[np.random.choice(range(len(sample_sizes)))] -= 1
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=k)
        for k in sample_sizes]
    # run the simulation
    begin_time = time.time()
    ts = msprime.simulate(
            population_configurations=population_configurations,
            migration_matrix=M,
            **kwargs)
    end_time = time.time()
    print("Simulation took {} seconds.".format(end_time - begin_time))
    return ts


def write_sim(filename, x, statname):
        outfl = open("{}.{}.tsv".format(filename, statname), "w")
        outfl.write("branch_{}\tsite_{}\n".format(statname, statname))
        for i in range(x.shape[0]):
            outfl.write("{}\t{}\n".format(x[i, 0], x[i, 1]))
        outfl.close()
        print("Written out to {}.{}.tsv".format(filename, statname))


if __name__ == "__main__":
    np.random.seed(111)
    n = 10
    # stepping stone migration
    # outmigration rate
    outmig = 0.0001
    M = np.zeros(shape=(n, n))
    for ii in range(n):
        M[ii, ii] = 0
        if ii > 0:
            M[ii, ii - 1] = 1 * outmig / 2
        if ii < n - 1:
            M[ii, ii + 1] = 1 * outmig / 2

    ts = tuple(run_sim(M,
                       nsamples=500,
                       length=1e6, Ne=1e4,
                       recombination_rate=1e-8,
                       mutation_rate=1e-8,
                       num_replicates=5))
    tsl = tuple(run_sim(M,
                        nsamples=500,
                        length=1e6, Ne=1e4,
                        recombination_rate=1e-8,
                        mutation_rate=1e-10,
                        num_replicates=5))
    nstats = int(n * (n-1) / 2 / 2)
    for statname in ("Y3", "f4"):
        print("Computing {} statistics: {}".format(nstats, statname))
        begin_time = time.time()
        x = np.concatenate(tuple(lmsp.compute_stat(t, statname, nstats) for t in ts))
        xl = np.concatenate(tuple(lmsp.compute_stat(t, statname, nstats) for t in tsl))
        end_time = time.time()
        print("Computing statistics took {} seconds.".format(end_time - begin_time))
        write_sim('simple', x, statname)
        write_sim('simple-lowmut', xl, statname)
