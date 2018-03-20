import msprime
import numpy as np
import landscape_msprime


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
            time=t, initial_size=x, growth_rate=0)
        for t, x in zip(change_times, change_factors)]

    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=M,
        demographic_events=demographic_events)
    dd.print_history()

    # run the simulation
    replicates = msprime.simulate(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=M,
        num_replicates=num_replicates)

    # And then iterate over these replicates
    T = np.zeros(num_replicates)
    for i, tree_sequence in enumerate(replicates):
        tree = tree_sequence.first()
        T[i] = tree.time(tree.root) / 4


if __name__ == "__main__":
    migration_example()
