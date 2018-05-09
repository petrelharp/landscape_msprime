import msprime
import numpy as np
import landscape_msprime

# TODO: what to do about empty cells? This is a hack:
MIN_SIZE = 100

def read_tree_nums():
    trees = np.loadtxt("tree_num_epsg3035.tsv")
    return trees

def read_coords():
    locfile = open("coords_epsg3035.tsv", "r")
    header = locfile.readline().split()
    assert(header == ['"X"', '"Y"', '"pop"', '"id"', '"cell"', '"msp_id"'])
    poplist = []
    idlist = []
    celllist = []
    for line in locfile:
        x, y, pop, id, cell, msp = [u.strip('"\n') for u in line.split()]
        poplist.append(pop)
        idlist.append(id)
        celllist.append(int(cell))
    return poplist, idlist, celllist


def population_expansion(num_replicates=1):
    tree_num = read_tree_nums()
    M = landscape_msprime.read_migration_matrix("abiesalba_epsg3035.migr.tsv")
    n = M.shape[0]
    ngens = tree_num.shape[0]
    assert(n == M.shape[1])
    assert(n == tree_num.shape[1])

    # obtain sample numbers from coordinates
    # recall we've ordered the coordinates so that msprime sample numbers
    # match the input
    poplist, idlist, celllist = read_coords()
    nsamples = [0 for _ in range(n)]
    for k in range(n):
        nsamples[k] = sum([u == k for u in celllist])

    population_configurations = [
        msprime.PopulationConfiguration(sample_size=nsamples[k],
            initial_size = max(MIN_SIZE, tree_num[ngens-1,k]))
        for k in range(n) ]

    # TODO: check the direction of time in tree_num ?!?!?
    # step forwards through history
    demographic_events = []
    t = ngens - 1
    x = tree_num[t,]
    while t > 0:
        t -= 1
        next_x = tree_num[t,]
        # the number of migrants from i to j is (M * x)[i,j]
        # and so the probability a lineage at j goes to i is (M*x)[i,j]/sum_k (M*x)[k,j]
        N = M * x
        N = N/N.sum(axis=0)
        N[np.isnan(N)] = 0.0
        for i in range(n):
            N[i,i] = 0.0
            if next_x[i] != x[i]:
                demographic_events.append(
                    msprime.PopulationParametersChange(
                        time=10*t, initial_size=max(MIN_SIZE, x[i]), 
                        population_id=i, growth_rate=0))
        demographic_events.append(
            msprime.MigrationRateChange(
                time=10*t, rate=N))
        x = next_x

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
