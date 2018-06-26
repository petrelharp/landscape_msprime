#!/usr/bin/env python3
description = '''
Simulate a coalescent within the history of Abies expansion.
'''

import msprime
import numpy as np
import os
import landscape_msprime

# for debugging: put pdb.set_trace() to drop into interactive
import pdb

### command-line parsing: execute this with '-h' for help
import argparse

parser = argparse.ArgumentParser(description=description)
parser.add_argument("--tree_nums", "-t", type=str, dest="tree_nums_file", required=True,
                    help="name of file to load tree number history from")
parser.add_argument("--sample_coords", "-s", type=str, dest="sample_coords_file", required=True,
                    help="name of file containing locations of sampling coordinates")
parser.add_argument("--migr_mat", "-m", type=str, dest="migr_mat_file", required=True,
                    help="name of file containing migration matrix")
parser.add_argument("--basedir", "-o", type=str, dest="basedir", required=True,
                    help="name of directory to save output files to.")
parser.add_argument("--dt", "-d", type=float, dest="dt", required=True,
                    help="time interval between steps")
parser.add_argument("--recomb", "-r", type=float, dest="recomb", required=True,
                    help="total recombination rate")
parser.add_argument("--mut", "-u", type=float, dest="mut", required=True,
                    help="total mutation rate")
parser.add_argument("--num_reps", "-n", type=float, dest="num_reps", default=1,
                    help="number of replicates")
parser.add_argument("--logfile", "-g", type=float, dest="logfile", 
                    help="name of log file")

args = parser.parse_args()

if not os.path.isdir(args.basedir):
    os.mkdir(args.basedir)

if args.logfile is None:
    args.logfile = os.path.join(args.basedir, "sim_trees.log")

logfile = open(args.logfile, "w")

# TODO: what to do about empty cells? This is a hack:
MIN_SIZE = 100

def read_tree_nums():
    # e.g. tree_num_epsg3035.tsv
    logfile.write("Reading tree numbers from {}\n".format(args.tree_nums_file))
    trees = np.loadtxt(args.tree_nums_file)
    return trees


def read_coords():
    # e.g. coords_epsg3035.tsv
    logfile.write("Reading sample coordinates from {}\n".format(args.sample_coords_file))
    locfile = open(args.sample_coords_file, "r")
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

def write_data(ts, key):
    # write out (a) the tree sequence, and (b) the haplotypes
    treefile = os.path.join(args.basedir, "rep.{}.trees".format(key))
    ts.dump(treefile)
    hapfile = os.path.join(args.basedir, "haps.{}.tsv".format(key))
    # haplotypes are in a matrix with one column per individual and one row per variant
    haps = ts.genotype_matrix()
    np.savetxt(hapfile, haps, fmt="%d")

def population_expansion(num_replicates=1):
    tree_num = read_tree_nums()
    # e.g. abiesalba_epsg3035.migr.tsv
    logfile.write("Reading migration matrix from {}\n".format(args.migr_mat_file))
    M = landscape_msprime.read_migration_matrix(args.migr_mat_file)
    n = M.shape[0]
    ngens = tree_num.shape[0]
    assert(n == M.shape[1])
    assert(n == tree_num.shape[1])
    logfile.write("Will run for {} generations, on {} populations.\n".format(ngens, n))

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
    last_N = M
    while t >= 0:
        t_ago = args.dt * (ngens - t)
        t -= 1
        next_x = tree_num[t,]
        # the number of migrants from i to j is (M * x)[i,j]
        # and so the probability a lineage at j goes to i is (M*x)[i,j]/sum_k (M*x)[k,j]
        N = (M.T * x).T
        N = N/N.sum(axis=0)
        N[np.isnan(N)] = 0.0
        for i in range(n):
            N[i,i] = 0.0
            if next_x[i] != x[i]:
                demographic_events.append(
                    msprime.PopulationParametersChange(
                        time=t_ago, initial_size=max(MIN_SIZE, x[i]), 
                        population_id=i, growth_rate=0))
        for i in range(n):
            for j in range(n):
                if (i != j) and (N[i,j] != last_N[i,j]):
                    demographic_events.append(
                        msprime.MigrationRateChange(
                            time=t_ago, rate=N[i,j], matrix_index=(i,j)))
        x = next_x
        last_N = N

    for x in demographic_events:
        logfile.write(str(x.time) + " :: " + x.__str__() + "\n")

    # run the simulation
    replicates = msprime.simulate(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=M,
        num_replicates=args.num_reps,
        recombination_rate=args.recomb,
        mutation_rate=args.mut)

    for j, rep in enumerate(replicates):
        write_data(rep, j)
        # # uncomment this to see all the genealogical trees:
        # for t in rep.trees():
        #     print(t.draw(format='unicode'))

    return True


if __name__ == "__main__":
    population_expansion()
