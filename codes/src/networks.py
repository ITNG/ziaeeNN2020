import matplotlib.pyplot as plt
import os
import igraph
import numpy as np
import pylab as pl
from sys import exit
from copy import copy
import networkx as nx
pl.switch_backend('agg')


dirName = ['../data/fig/', '../data/text/', '../data/npz/', 'networks']
for d in dirName:
    if not os.path.exists(d):
        os.mkdir(d)

# np.set_printoptions(threshold=np.nan)


class make_graph:
    ''' make different graphs ans return their adjacency matrices
    as a 1 dimensional double vector in stl library'''

    def __init__(self, seed=None):

        self.G = 0
        if seed:
            self.seed = seed
            np.random.seed(seed)

    # ---------------------------------------------------------------#
    def print_adjacency_matrix(self, G):
        '''print adjacency matrix of graph on the screen '''
        M = nx.to_numpy_matrix(G)
        M = np.array(M)
        for row in M:
            for val in row:
                print ('%.0f' % val,)
            print("")

    # ---------------------------------------------------------------#
    def complete_graph(self, N, plot_adj=False):
        ''' returns all to all graph '''

        self.N = N
        self.G = nx.complete_graph(N)
        M = nx.to_numpy_matrix(self.G)

        if plot_adj:
            self.imshow_plot(M, "con")

        return np.asarray(M)
    # ---------------------------------------------------------------#

    def complete_weighted(self, N, clusters, weights, d, plot_adj=False):
        '''
        return a complete weighted graph, weights distribute in clusters
        with the same size.
        weights: list with length 2, shows the weights of edges
        in and between clusters, respectivly
        clusters : list of cluster lengths
        d : delay for all nodes
        '''

        self.N = N
        # I use delay as weight
        self.modular_graph(N, 1, 1, clusters, weights[0], weights[1])
        M = self.D
        if plot_adj:
            self.imshow_plot(np.asarray(M).reshape(N, N), "con")

        self.D = self.complete_graph(N) * d
        return M
    # ---------------------------------------------------------------#

    def complete_hmn_weighted(self, n_M0, level, weights, delay,
                              plot_adj=False, seed=124):
        '''
        return a complete weighted graph, weights distribute in a hierarchical
        modular form.
        Single delay for every node.
        '''
        self.hierarchical_modular_graph(n_M0, level, 1, 1, 1, weights,
                                        plot_adj=False, seed=124)
        N = self.N
        M = nx.to_numpy_matrix(self.G, weight='weight')

        if plot_adj:
            self.imshow_plot(M, "con")
        self.D = (delay,)*N*N

        return tuple(np.asarray(M).reshape(-1))

    # ---------------------------------------------------------------#
    def erdos_renyi_graph(self, N, p,
                          seed=123,
                          directed=False,
                          plot_adj=False):
        ''' returns Erdos Renyi graph '''

        self.N = N
        self.G = nx.erdos_renyi_graph(N, p, seed=seed, directed=False)
        M = nx.to_numpy_matrix(self.G)

        if plot_adj:
            self.imshow_plot(M, "con")

        return M
    # ---------------------------------------------------------------#

    def modular_graph(self, N, pIn, pOut, lengths, dIn, dOut, plot_adj=False):
        ''' returns a modular networks with :
        N : number of nodes
        pIn : conectivity of nodes insede clusters
        pOut: conectivity of nodes between clusters
        n_cluster :  number of clusters in graph
        dIn  : delay between nodes inside the clusters
        dOut : delay between nodes outside the clusters
        '''

        self.N = N
        M = np.zeros((N, N))
        D = np.zeros((N, N))
        n_cluster = len(lengths)

        for i in range(N):
            for j in range(i+1, N):
                r = np.random.rand()
                if r < pOut:
                    M[i, j] = M[j, i] = 1.0
                    D[i, j] = D[j, i] = dOut

        # empty inside the clusters
        s = 0
        for k in range(n_cluster):
            if k > 0:
                s += lengths[k-1]
            for i in range(s, (lengths[k]+s)):
                for j in range(i+1, (lengths[k]+s)):
                    M[i, j] = M[j, i] = 0.0
                    D[i, j] = D[j, i] = 0.0

        # fill inside the clusters
        s = 0
        for k in range(n_cluster):
            if k > 0:
                s += lengths[k-1]
            for i in range(s, (lengths[k]+s)):
                for j in range(i+1, (lengths[k]+s)):
                    r = np.random.rand()
                    if r < pIn[k]:
                        M[i, j] = M[j, i] = 1.0
                        D[i, j] = D[j, i] = dIn

        # print delay matrix
        def print_delay_matrix():
            ofi = open("delay.txt", "w")
            for i in range(N):
                for j in range(N):
                    ofi.write("%2.0f" % D[i, j])
                ofi.write("\n")
            ofi.close()

        self.G = nx.from_numpy_matrix(M)
        self.D = D

        if plot_adj:
            self.imshow_plot(M, "con")
            # self.print_adjacency_matrix(self.G)

        return M
    # ---------------------------------------------------------------#

    def hierarchical_modular_graph(self, n_M0, level, prob0, prob, alpha,
                                   delays, plot_adj=False):
        '''
        n_M0 : size of module at level 1
        s    : number of levels
        n_modules : number of modules
        N : number of nodes
        ps : probability of conection in each level
             level one is 1 and the others determine with prob function
        delays: delay in each level as a list
        '''
        def probability(l, a=1, p=0.25):
            if l == 0:
                from sys import exit
                print ("level shold be integer and > 0", l)
                exit(0)
            else:
                return a * p**l

        s = level
        n_modules = int(2**(s-1))  # number of modules
        N = int(n_modules*n_M0)  # number of nodes
        self.N = N

        # M0 = nx.complete_graph(n_M0)
        M0 = nx.erdos_renyi_graph(n_M0, prob0, seed=self.seed)
        for e in nx.edges(M0):
            M0.add_edge(*e, weight=delays[0])  # delays act in weight attribute

        ps = [prob0]+[probability(i, alpha, prob) for i in range(1, s)]

        for l in range(1, s):
            if l == 1:
                M_pre = [M0] * n_modules
            else:
                M_pre = copy(M_next)
            M_next = []
            k = 0
            for ii in range(n_modules/(2**l)):
                step = 2**(l-1)
                tG = nx.convert_node_labels_to_integers(M_pre[k+1], step*n_M0)
                tG1 = nx.compose(M_pre[k], tG)
                edge = 0
                effort = 1
                ''' make sure that connected modules are not isolated '''
                while edge < 1:
                    # print "effort ", effort
                    effort += 1
                    for i in range(len(tG1)):
                        for j in range(i+1, len(tG1)):
                            if (i < step*n_M0) & (j > step*n_M0-1) & (np.random.rand() < ps[l]):
                                tG1.add_edge(i, j, weight=delays[l])
                                edge += 1

                M_next.append(tG1)
                k += 2
        self.G = M_next[0]

        if plot_adj:
            M = nx.to_numpy_matrix(self.G, weight=None)
            self.imshow_plot(M, "con")
            # self.print_adjacency_matrix(self.G)
        D = nx.to_numpy_matrix(self.G, weight='weight')
        self.D = np.asarray(D)

        M = nx.to_numpy_matrix(self.G, weight=None)

        return np.asarray(M)
    # ---------------------------------------------------------------#

    def from_adjacency_matrix_graph(self, filename, plot_adj=False):
        ''' makes a graph from adjacency matrix in filename
        and return 1D double vector in stl '''

        A = np.genfromtxt(filename)
        self.N = len(A)

        if plot_adj:
            self.imshow_plot(A, "con")

        return A
    # ---------------------------------------------------------------#

    def imshow_plot(self, data, name, title=None):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        fig = pl.figure(140, figsize=(6, 6))
        pl.clf()
        ax = pl.subplot(111)
        im = ax.imshow(data, interpolation='nearest',
                       cmap='afmhot')  # , cmap=pl.cm.ocean
        # ax.invert_yaxis()
        if title:
            ax.set_title(title)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        pl.colorbar(im, cax=cax)
        pl.savefig('../data/fig/'+name+".png", dpi=300)
        # pl.close()
    # ---------------------------------------------------------------#

    def multilevel(self, data):
        conn_indices = np.where(data)
        weights = data[conn_indices]
        edges = zip(*conn_indices)
        G = igraph.Graph(edges=edges, directed=False)
        G.es['weight'] = weights
        comm = G.community_multilevel(weights=weights, return_levels=False)
        return comm
    # ---------------------------------------------------------------#

    def hmn_cortex_alike(self, len_communities=None,
                         pIn_clusters=None, p2=0.2, p3=0.04,
                         weight1=1.0, weight2=0.5, weight3=0.2,
                         delay1=5.0, delay2=20.0, delay3=25.0,
                         plot_adj=False):

        if len_communities is None:
            len_communities = [12, 13, 14, 14, 12]
            pIn_clusters = [0.515, 0.731, 0.670, 0.527, 0.712]

        number_of_communities = len(len_communities)
        assert (number_of_communities == len(pIn_clusters))
        N = np.sum(len_communities)

        # ---------------------------------------------------------------#
        def fill_in_modules():

            G = nx.Graph()
            M = []
            for i in range(number_of_communities):
                tg = nx.erdos_renyi_graph(len_communities[i],
                                          pIn_clusters[i], seed=self.seed)
                M.append(tg)
                M[i] = nx.convert_node_labels_to_integers(M[i], nodes[i][0])
                G = nx.compose(G, M[i])
            return G.edges()
        # ---------------------------------------------------------------#

        def fill_between_modules(p2, p3):

            def empty_inside_modules(M):
                s = 0
                for k in range(number_of_communities):
                    if k > 0:
                        s += len_communities[k-1]
                    for i in range(s, (len_communities[k]+s)):
                        for j in range(i+1, (len_communities[k]+s)):
                            M[i, j] = M[j, i] = 0.0
                return M

            M = np.zeros((N, N))

            for i in range(N/2):
                for j in range(i+1, N/2):
                    r = np.random.rand()
                    if r < p2:
                        M[i, j] = M[j, i] = 1.0

            for i in range(N/2, N):
                for j in range(i+1, N):
                    r = np.random.rand()
                    if r < p2:
                        M[i, j] = M[j, i] = 1.0

            M = empty_inside_modules(M)
            tg2 = nx.from_numpy_matrix(M)
            edges2 = tg2.edges()

            M = np.zeros((N, N))
            for i in range(N):
                for j in range(i+1, N):
                    r = np.random.rand()
                    if r < p3:
                        M[i, j] = M[j, i] = 1.0

            M = empty_inside_modules(M)
            tg3 = nx.from_numpy_matrix(M)
            edges3 = tg3.edges()

            return edges2, edges3
        # ---------------------------------------------------------------#
        # pOut = 0.157

        nodes = []
        a = [0]
        for i in range(number_of_communities):
            a.append(a[i]+len_communities[i])
            nodes.append(range(a[i], a[i+1]))

        edges1 = fill_in_modules()
        edges2, edges3 = fill_between_modules(p2, p3)

        G = nx.Graph()
        for e in edges1:
            G.add_edge(*e, weight=weight1, delay=delay1)
        for e in edges2:
            G.add_edge(*e, weight=weight2, delay=delay2)
        for e in edges3:
            G.add_edge(*e, weight=weight3, delay=delay3)

        if plot_adj:
            A = extract_attributes(G, "weight")
            D = extract_attributes(G, "delay")
            self.imshow_plot(A, "A", title="Coupling weights")
            self.imshow_plot(D, "D", title='Delays')

        print (nx.info(G))

        return A, D

        # ---------------------------------------------------------------#


def extract_attributes(G, attr=None):
    edges = G.edges()
    n = nx.number_of_nodes(G)
    A = np.zeros((n, n))

    for e in edges:
        A[e[0], e[1]] = A[e[1], e[0]] = G[e[0]][e[1]][attr]
    return A
    # --------------------------------------------------------------#


if __name__ == "__main__":
    np.random.seed(1235)
    import os
    clusters = [20, 20, 20]
    gr = make_graph()
    Adj = gr.modular_graph(sum(clusters),
                           0.7, 0.1, clusters,
                           2.0, 4.3, True)
    Dij = gr.D
    np.savetxt("dat/C.dat", Adj, fmt="%d")
    np.savetxt("dat/D.dat", Dij, fmt="%15.6f")
