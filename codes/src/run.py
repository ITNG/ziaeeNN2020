import os
import lib
import numpy as np
import networks
from os import system
from time import time
from threading import Thread
from joblib import Parallel, delayed


# preparing directories
directories = ['../data/fig', '../data/npz', '../data/text']
for d in directories:
    if not os.path.exists(d):
        os.makedirs(d)


def run_command(arg):
    
    command = "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(*arg)
    system("./prog " + command)


def batch_run():
    arg = []
    for g in G:
        for m in mu:
            arg.append([N,
                        t_trans,
                        t_sim,
                        g,
                        m,
                        cname,
                        dname,
                        ncluster,
                        num_sim,
                        cor_to_file])
    Parallel(n_jobs=n_jobs)(
        map(delayed(run_command), arg))

# parameters--------------------------------------------------------#

clusters1 = [12, 13, 14, 14, 12]        # size of clusters in level 1
clusters2 = [32, 33]                    # size of clusters in level 2
ncluster = len(clusters1)               # number of clusters in level 1
N = sum(clusters1)                      # number of nodes
G = [3.0]                               # coupling (it will divided by G/N in the simulations)
mu = np.arange(0.0, 1.5, 0.01)          # inital frequencies
num_sim = 50                            # number of ensembles
t_trans = 3000.0                        # transition time in [ms]
t_sim = t_trans + 7005.0                # simulation time in [ms]
cname = "C"                             
dname = "D"
cor_to_file = 1                         # flag to record the correlations in file
n_jobs = 4                              # number of process for parallel simulations (number of avilable cores)
# ------------------------------------------------------------------#

if __name__ == "__main__":

    start = time()
    seed = 136

    gr = networks.make_graph(seed)
    Adj = np.loadtxt("networks/C65.txt")
    Dij = np.loadtxt("networks/L65.txt")/5.0
    Adj = Adj/np.max(Adj)

    np.savetxt("networks/C.txt", Adj, fmt="%15.6f")
    np.savetxt("networks/D.txt", Dij, fmt="%15.6f")

    batch_run()
    lib.display_time(time()-start)