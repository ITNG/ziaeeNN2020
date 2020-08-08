import numpy as np
import networkx as nx
import igraph
import pylab as pl
from scipy.stats.stats import pearsonr


adrsf = "../data/f/"
adrs = "../data/"

# -------------------------------------------------------------------#


def walktrap(adj, steps=4):
    conn_indices = np.where(adj)
    weights = adj[conn_indices]
    edges = list(zip(*conn_indices))
    G = igraph.Graph(edges=edges, directed=False)
    comm = G.community_walktrap(weights, steps=steps)
    communities = comm.as_clustering()
    # print comm
    # print "clusters : ", communities
    # print "optimal count : ", comm.optimal_count
    return communities
# -------------------------------------------------------------------#


def comm_weighted(data, remove_self_loops=True):
    '''Community structure based on the multilevel
       algorithm of Blondel et al.'''

    conn_indices = np.where(data)
    weights = data[conn_indices]
    edges = list(zip(*conn_indices))
    if remove_self_loops:
        for e in edges:
            if e[0] == e[1]:
                edges.remove(e)
    G = igraph.Graph(edges=edges, directed=False)
    G.es['weight'] = weights
    comm = G.community_multilevel(weights=weights, return_levels=False)
    return comm
# -------------------------------------------------------------------#


def comm_unweighted(data, remove_self_loops=True):
    '''Community structure based on the multilevel
       algorithm of Blondel et al.'''

    conn_indices = np.where(data)
    edges = zip(*conn_indices)
    if remove_self_loops:
        for e in edges:
            if e[0] == e[1]:
                edges.remove(e)
    G = igraph.Graph(edges=edges, directed=False)
    comm = G.community_multilevel(weights=None, return_levels=False)
    return comm
# -------------------------------------------------------------------#


def calculate_NMI(comm1, comm2):
    '''Compares two community structures using normalized
    mutual information as defined by Danon et al (2005)'''

    nmi = igraph.compare_communities(
        comm1, comm2, method='nmi', remove_none=False)
    return nmi

# -------------------------------------------------------------------#


def display_time(time):
    ''' '''

    hour = int(time/3600)
    minute = (int(time % 3600))//60
    second = time-(3600.*hour+60.*minute)
    print ("Done in %d hours %d minutes %09.6f seconds" \
        % (hour, minute, second))
# -------------------------------------------------------------------#


def binarize(data, threshold):
    data = np.asarray(data)
    upper, lower = 1, 0
    data = np.where(data >= threshold, upper, lower)
    return data
# -------------------------------------------------------------------#


def imshow_plot(data, fname='R', cmap='afmhot',
                figsize=(5, 5),
                vmax=None, vmin=None):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    fig = pl.figure(100, figsize=figsize)
    pl.clf()
    ax = pl.subplot(111)
    im = ax.imshow(data, interpolation='nearest', cmap=cmap,
                   vmax=vmax, vmin=vmin)
    ax.invert_yaxis()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pl.colorbar(im, cax=cax)
    pl.savefig(fname, dpi=150)
    pl.close()

# ---------------------------------------------------------------------- #


def reorder_nodes(C, communities):

    # reordering the nodes:
    N = C.shape[0]
    n_comm = len(communities)

    nc = []
    for i in range(n_comm):
        nc.append(len(communities[i]))

    # new indices of nodes-------------------------------------------

    newindices = []
    for i in range(n_comm):
        newindices += communities[i]
    # --------------------------------------------------------------

    reordered = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            reordered[i, j] = C[newindices[i], newindices[j]]

    return reordered
# ---------------------------------------------------------------------- #


def read_from_file(filename):
    """
    Generate lists using data from file.

    :param filename:
    :return: list of data (list of list)
    """
    DATA = []
    with open(filename, "r") as datafile:
        lines = datafile.readlines()
        for line in lines:
            data = line.split()
            data = [int(i) for i in data]
            DATA.append(data)

    return DATA
# ---------------------------------------------------------------------- #


def nodes_of_each_cluster(len_communities):
    "return nodes of each cluster as list of list"

    n_comm = len(len_communities)
    nodes = []
    a = [0]
    for i in range(n_comm):
        a.append(a[i]+len_communities[i])
        nodes.append(range(a[i], a[i+1]))
    return nodes
