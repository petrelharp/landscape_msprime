import msprime
import numpy as np
import time

import landscape_msprime as lmsp

def run_sim(migr_matrix_file, nsamples, **kwargs):
    M = lmsp.read_migration_matrix(migr_matrix_file)
    assert(M.shape[0] == M.shape[1])
    n = M.shape[0]
    print("Simulating on a landscape of {} patches.".format(n))
    # sample a total of nsamples, uniformly spread
    sample_sizes = [int(np.ceil(nsamples/n)) for _ in range(n)]
    while sum(sample_sizes) < nsamples:
        sample_sizes[np.random.choice(range(len(samples)))] -= 1
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



if __name__ == "__main__":
    ts = run_sim("nussear_8x.migr.tsv", 
            nsamples=500, 
            length=1e6, Ne=1e4,
            recombination_rate=1e-8,
            mutation_rate=1e-8)
    nstats = 100
    for statname in ("f2", "f3", "f4", "Y2", "Y3"):
        print("Computing {} statistics: {}".format(nstats, statname))
        begin_time = time.time()
        x = lmsp.compute_stat(ts, statname, nstats)
        end_time = time.time()
        print("Computing statistics took {} seconds.".format(end_time - begin_time))
        outf = open("nussear_8x.{}.tsv".format(statname), "w")
        outf.write("branch_{}\tsite_{}\n".format(statname, statname))
        for i in range(x.shape[0]):
            outf.write("{}\t{}\n".format(x[i,0], x[i,1]))
        outf.close()
        print("Written out to nussear_8x.{}.tsv".format(statname))
