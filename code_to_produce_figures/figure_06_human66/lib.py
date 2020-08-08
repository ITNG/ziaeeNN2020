import numpy as np
import pylab as pl
import networkx as nx
from mpl_toolkits.axes_grid1 import make_axes_locatable


fontsize = 14
labelsize = 14


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
        ax.set_xlabel(r"$\sigma$", fontsize=16, labelpad=0)
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
    if label:
        label = label / (2.0 * np.pi)*1000
    axs.text(0.05, 0.3, str(r"$\nu_0 = %.0f$Hz" % label),
             fontsize=14, transform=axs.transAxes)
    axs.tick_params(labelsize=14)
    axs.text(-0.3, 0.93, label1,
             fontsize=14, transform=axs.transAxes)


def plot_rectangles(comm, ordered, ax,
                    title=None, ylabel=False):

    from matplotlib.patches import Rectangle

    if title:
        ax.set_title(title, fontsize=14)
    s = 0
    X = Y = 0
    N = ordered.shape[0]
    n_comm = len(comm)
    for k in range(n_comm):
        if k > 0:
            s += len(comm[k-1])
            X = s
            Y = s
        ax.add_patch(Rectangle((X-0.5, Y-0.5),
                               len(comm[k]), len(comm[k]),
                               fill=None, lw=1.5, alpha=0.5,))
        # alpha=1,linestyle='--',lw=3,facecolor='none'))

    im = ax.imshow(ordered, interpolation='nearest',
                   cmap="jet")
    ax.invert_yaxis()
    cax = ax.inset_axes([1, 0.0, 0.05, 0.4])
    cbar = pl.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=labelsize-2)
    ax.axvline(x=32.5, ls="--", lw=2, c='k', alpha=0.8)
    ax.axhline(y=32.5, ls="--", lw=2, c='k', alpha=0.8)

    ax.set_xticks([])
    ax.set_yticks([])
    if ylabel:
        ax.set_ylabel("region", fontsize=14)
    ax.set_xlabel("region", fontsize=14)
# ------------------------------------------------------- #


def read_from_file(filename, n):
    Data = []
    with open(filename, "r") as f:
        for i in range(n):
            data = f.readline().split()
            data = [int(i) for i in data]
            Data.append(data)
    return Data
# ------------------------------------------------------- #


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx  # , array[idx]
# ------------------------------------------------------- #


def find_first_value(arr, value, tol=0.01):
    for i, item in enumerate(arr):
        if np.abs(item-value) < tol:
            return i
    print ("match not found!")
    exit(0)
# ------------------------------------------------------- #


def toRad(x):
    return x*(2.*np.pi)/1000.0
# ------------------------------------------------------- #


def toHz(x):
    return x/(2.0*np.pi)*1000
# ------------------------------------------------------- #


def cor_plot(ax, cor, vmax=None, vmin=None,
             lb=None, ylabel=False,
             ticks=None, format='%.0f'):

    im = ax.imshow(cor,
                   interpolation='nearest',
                   cmap='afmhot',
                   alpha=0.9,
                   vmin=vmin, vmax=vmax,
                   aspect="auto"
                   )
    ax.invert_yaxis()
    divider = make_axes_locatable(ax)
    cax = ax.inset_axes([1, 0.0, 0.05, 0.3])
    cbar = pl.colorbar(im, cax=cax, ticks=ticks)
    cbar.ax.tick_params(labelsize=labelsize-2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(-0.25, 0.85, lb,
            fontsize=fontsize,
            transform=ax.transAxes)
    ax.set_xlabel("index", fontsize=14)
    if ylabel:
        ax.set_ylabel("index", fontsize=14)
# ------------------------------------------------------- #
