import netCDF4
import msprime
import numpy as np
import landscape_msprime


def read_tree_nums():
    trees = netCDF4.Dataset("abiesalba_h1x1.cdf", "r", format="NETCDF4")
    tree_num = np.zeros(shape=(2204, 55, 70))
    for t in range(2204):
        # recall python is zero-indexed so '9' = 10th variable in fpc_grid
        tree_num[t,:,:] = trees['lu_area'][t,0,:,:] * trees['area'] * trees['fpc_grid'][t,9,:,:]/40.
    # the result has NAs as -99999
    tree_num[tree_num < 0] = 0.0
    return tree_num


def population_expansion(num_replicates=1):
    tree_num = read_tree_nums()
    M = landscape_msprime.read_migration_matrix("abiesalba_longlat.migr.tsv")
    assert(M.shape[0] == M.shape[1])
    n = M.shape[0]

    # sample one individual per deme
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=1)
        for _ in range(n) ]
    demographic_events = []

    # step forwards through history
    for t in range(tree_num.shape[2],0,-1):
        x = tree_num[,,t]
        # the number of migrants from i to j is (M * x)[i,j]
        # and so the probability a lineage at j goes to i is (M*x)[i,j]/sum_k (M*x)[k,j]
        N = M * x.flatten()
        N = N/N.sum(axis=0)
        N[np.isnan(N)] = 0.0
        for i in range(n):
            N[i,i] = 0.0

        demographic_events.append(
            msprime.PopulationParametersChange(
                time=10*t, initial_size=x, growth_rate=0))
        demographic_events.append(
            msprime.MigrationRateChange(
                time=10*t, rate=N))

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
    population_expansion()
