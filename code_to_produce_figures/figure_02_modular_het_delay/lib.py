import numpy as np
import pylab as pl
import networkx as nx


def nodes_of_each_cluster(len_communities):
    "return nodes of each cluster as list of list"

    n_comm = len(len_communities)
    nodes = []
    a = [0]
    for i in range(n_comm):
        a.append(a[i]+len_communities[i])
        nodes.append(range(a[i], a[i+1]))
    return nodes
# ----------------------------------------- #


def divide_edges_in_between_clusters(
        mat,
        communities):
    """
    divide edges in and between clusters.
    print two matrix.
    """
    n_comm = len(communities)
    N = mat.shape[0]

    len_communities = []
    for i in range(n_comm):
        len_communities.append(len(communities[i]))

    nodes = nodes_of_each_cluster(len_communities)
    G = nx.from_numpy_matrix(mat)
    edges = G.edges()
    list_in = []
    list_between = []

    for e in edges:
        if e[0] == e[1]:
            continue
        info = check_edge_between_clusters(e, nodes)
        if info:
            list_between.append(mat[e[0], e[1]])
        else:
            list_in.append(mat[e[0], e[1]])

    return list_in, list_between


def check_edge_between_clusters(e, nodes):
    "check if given edge is between clusters."
    # print e[0], e[1]
    n_clusters = len(nodes)
    for i in range(n_clusters):
        if ((e[0] in nodes[i]) & (e[1] not in nodes[i])) |\
                ((e[0] not in nodes[i]) & (e[1] in nodes[i])):
            return True

    return False


def plot_hist(arr, num_of_points, ax, bins='auto', color='red',
              label=None,
              xticks=None,
              yticks=None,
              xlim=None,
              ylim=None,
              xlabel=None,
              ylabel=None,
              logscale=False,
              bin_range=(-1, 1)):

    hist, bins = np.histogram(arr, bins=bins, range=bin_range)
    hist = hist / float(num_of_points)

    width = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.bar(center, hist, align='center',
           width=width,
           alpha=0.7,
           color=color,
           label=label)
    ax.legend(frameon=False, loc='upper left', fontsize=14)
    if logscale:
        ax.set_yscale('log')
    if xticks:
        ax.set_xticks(xticks)
    if yticks:
        ax.set_yticks(yticks)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if xlabel:
        ax.set_xlabel(r"$\sigma$", fontsize=16)
    if ylabel:
        ax.set_ylabel("Prob", fontsize=14, labelpad=0)
# ----------------------------------------- #


def plot_correlation_distribution(
        cor,
        clusters,
        axs,
        nbins=100,
        label=None,
        label1=None,
        thr=1e-10,
        ** kwargs):

    communities = nodes_of_each_cluster(clusters)
    corIn, corOut = divide_edges_in_between_clusters(cor, communities)
    num_of_points = (np.abs(cor) > thr).sum()

    plot_hist(corOut, num_of_points,  axs, bins=nbins, color='royalblue',
              label="inter", **kwargs)
    plot_hist(corIn, num_of_points, axs, bins=nbins, color='red',
              label="intra", **kwargs)
    axs.text(0.05, 0.4, str(r"$\nu_0 = %.0f$Hz" % label),
             fontsize=14, transform=axs.transAxes)
    axs.tick_params(labelsize=14)
    axs.text(-0.3, 0.91, label1,
             fontsize=14, transform=axs.transAxes)
