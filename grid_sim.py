import msprime
import numpy as np

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


def migration_example(num_replicates=1):
    M = read_migration_matrix("nussear_8x.migr.tsv")
    # M = read_migration_matrix("nussear.migr.tsv")
    assert(M.shape[0] == M.shape[1])
    n = M.shape[0]

    # sample one individual per deme
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=1)
        for _ in range(n) ]

    # suppose pop size changes have been sinusoidal
    change_times = range(100)
    change_factors = [1.0 + 0.5 * np.sin(2*np.pi*t/20) for t in change_times]
    demographic_events = [
        msprime.PopulationParametersChange(
            time=t, initial_size=change_factors, growth_rate=0, population_id=None),
        for t in change_times]

    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    dd.print_history()

    # run the simulation
    replicates = msprime.simulate(
        population_configurations=population_configurations,
        migration_matrix=M,
        num_replicates=num_replicates)

    # And then iterate over these replicates
    T = np.zeros(num_replicates)
    for i, tree_sequence in enumerate(replicates):
        tree = tree_sequence.first()
        T[i] = tree.time(tree.root) / 4


if __name__ == "__main__":
    migration_example()
