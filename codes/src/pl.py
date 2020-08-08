import os
import lib
import pylab as pl
import numpy as np
import pandas as pd
from sys import exit
from time import time
from run import N, G, mu, num_sim, clusters1, clusters2

# ---------------------------------------------------------------------- #


def imshow_plot(data):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    fig = pl.figure(100, figsize=(10, 10))
    pl.clf()
    ax = pl.subplot(111)
    im = ax.imshow(data, interpolation='nearest',
                   cmap='afmhot')  # , cmap=pl.cm.ocean
    ax.invert_yaxis()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pl.colorbar(im, cax=cax)

# ---------------------------------------------------------------------- #


def plot_phase_space(R, X, Y, name="R", xtickstep=1, ytickstep=1,
                     xlabel=None, ylabel=None, title=None,
                     vmax=None, vmin=None):
    '''
    plot R in 2D plane of X and Y axises
    '''
    print (len(X), len(Y), R.shape)

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    r, c = R.shape
    assert((r > 1) & (c > 1))

    x_step = X[1] - X[0]
    y_step = Y[1] - Y[0]

    f, ax = pl.subplots(1, figsize=(6, 6))
    im = ax.imshow(R, interpolation='nearest', aspect='auto',
                   cmap='afmhot', vmax=vmax, vmin=vmin)
    ax.invert_yaxis()

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pl.colorbar(im, cax=cax)

    ax.set_xticks(np.arange(0, len(X), xtickstep))
    ax.set_xticklabels(str("%.1f" % i)for i in X[::xtickstep])
    ax.set_yticks(np.arange(0, len(Y), ytickstep))
    ax.set_yticklabels(str("%.1f" % i)for i in Y[::ytickstep])

    if xlabel:
        ax.set_xlabel(xlabel, fontsize=16)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=16)
    if title:
        ax.set_title(title, fontsize=16)

    pl.savefig("fig/"+name+".png")
    pl.close()
# ---------------------------------------------------------------------- #


def plot_nmi(adj, G, Omega, threshold=0.95, to_npz=False):

    n = len(Omega)
    N = adj.shape[0]
    communities_l1 = lib.read_from_file(
        str('../src/networks/communities_%d_l1.txt' % N))
    communities_l2 = lib.nodes_of_each_cluster([32, 33])

    adj = lib.reorder_nodes(adj, communities_l1)

    adj_communities1 = lib.comm_weighted(adj)
    adj_membership1 = adj_communities1.membership
    adj_membership2 = [0]*32+[1]*33

    fig, ax = pl.subplots(1, figsize=(6, 5))
    for g in G:
        nmi1 = np.zeros((num_sim, n))
        nmi2 = np.zeros((num_sim, n))
        for i in range(n):
            for ens in range(num_sim):
                subname = str("%.6f-%.6f" % (g, Omega[i]))
                c = np.fromfile("text/c-"+subname+"-"+str(ens)+".bin",
                                dtype=float, count=-1)
                c = np.reshape(c, (N, N))
                c = lib.binarize(c, threshold)
                cor_communities1 = lib.walktrap(  # comm_unweighted
                    lib.reorder_nodes(c, communities_l1), steps=8)
                cor_communities2 = lib.walktrap(
                    lib.reorder_nodes(c, communities_l2), steps=60)
                # nmi[ens, i] = lib.calculate_NMI(adj_communities,
                #                                 cor_communities)
                nmi1[ens, i] = lib.calculate_NMI(
                    adj_membership1,
                    cor_communities1.membership)
                nmi2[ens, i] = lib.calculate_NMI(
                    adj_membership2,
                    cor_communities2.membership)
        if num_sim > 1:
            nmi_mean1 = np.mean(nmi1, axis=0)
            nmi_std1 = np.std(nmi1, axis=0)
            nmi_mean2 = np.mean(nmi2, axis=0)
            nmi_std2 = np.std(nmi2, axis=0)
        else:
            nmi_mean1 = nmi1[0, :]
            nmi_std1 = 0
            nmi_mean2 = nmi2[0, :]
            nmi_std2 = 0
        if to_npz:
            np.savez("npz/nmi-"+str('%.6f' % g),
                     nmi=[nmi_mean1, nmi_mean2],
                     std=[nmi_std1, nmi_std2],
                     omega=Omega,
                     g=g)

        if num_sim > 1:
            ax.errorbar(Omega, nmi_mean1, yerr=nmi_std1, lw=2,
                        marker='o', label=str("%.2f" % g))
            ax.errorbar(Omega, nmi_mean2, yerr=nmi_std2, lw=2,
                        marker='o', label=str("%.2f" % g))
        else:
            ax.plot(Omega, nmi_mean1, lw=2, label=str("%.2f" % g))
            ax.plot(Omega, nmi_mean2, lw=2, label=str("%.2f" % g))

    ax.legend()
    ax.set_xlabel(r"$\omega$", fontsize=14)
    ax.set_ylabel(r"NMI")
    pl.tight_layout()
    pl.savefig("fig/nmiq.png", dpi=300)
    pl.close()
# ---------------------------------------------------------------------- #


def plot_correlations(g, omega, ns=num_sim):

    n = len(omega)
    communities_l1 = lib.read_from_file(
        str('../src/networks/communities_%d_l1.txt' % N))
    for g in G:
        for i in range(n):
            cor = np.zeros((N, N))
            for ens in range(ns):
                subname = str("%.6f-%.6f" % (g, omega[i]))

                c = np.fromfile("text/c-"+subname+"-"+str(ens)+".bin",
                                dtype=float, count=-1)
                c = np.reshape(c, (N, N))
                cor += c
            cor /= float(ns)
            cor_ordered = lib.reorder_nodes(cor, communities_l1)
            np.savez("npz/cor-ave-"+subname, cor=cor_ordered)
            lib.imshow_plot(cor_ordered,
                            fname="fig/c-"+subname+'.png',
                            vmax=1, vmin=-1)
# ---------------------------------------------------------------------- #


def plot_r_local_global_hmn(G, mu, clusters1, clusters2):

    nx = len(mu)
    ny = len(G)
    clusters1 = np.asarray(clusters1)
    clusters2 = np.asarray(clusters2)
    prob1 = clusters1/float(np.sum(clusters1))
    prob2 = clusters2/float(np.sum(clusters2))
    nc = len(clusters1)

    Rg = np.zeros((ny, nx))
    Rl1 = np.zeros((ny, nx))
    Rl2 = np.zeros((ny, nx))

    for i in range(ny):
        for j in range(nx):

            fname = 'R-'+str('%.6f' % G[i])+'-'+str('%.6f' % mu[j])+'.txt'
            c = np.loadtxt('text/'+fname)
            if num_sim > 1:
                c = np.mean(c, axis=0)

            Rg[i, j] = c[2]
            rl1 = c[3:(3+nc)]
            rl2 = c[(3+nc):]
            Rl1[i, j] = np.sum(prob1*rl1)
            Rl2[i, j] = np.sum(prob2*rl2)

    np.savez('npz/R', Rg=Rg, Rl1=Rl1, Rl2=Rl2, g=G, mu=mu)

    if len(G) > 1:
        plot_phase_space(Rg, mu, G, 'Rg', xlabel='mu', ylabel='g',
                         xtickstep=5, ytickstep=10)
        plot_phase_space(Rl1, mu, G, 'Rl1', xlabel='mu', ylabel='g',
                         xtickstep=5, ytickstep=10)
        plot_phase_space(Rl2, mu, G, 'Rl2', xlabel='mu', ylabel='g',
                         xtickstep=5, ytickstep=10)
    else:
        fig, ax = pl.subplots(1, figsize=(8, 5))
        ax.plot(mu, Rg[0, :], lw=2, c='r', label=r'$R_g$')
        ax.plot(mu, Rl1[0, :], lw=2, c='b', label=r'$R_{l_1}$')
        ax.plot(mu, Rl2[0, :], lw=2, c='g', label=r'$R_{l_2}$')
        ax.legend()
        pl.savefig('fig/'+str('%.2f.png' % G[0]))

# ---------------------------------------------------------------------- #


def plot_r_local_global_modular(G, mu, clusters1):

    nx = len(mu)
    ny = len(G)
    clusters1 = np.asarray(clusters1)
    prob = clusters1/float(np.sum(clusters1))
    nc = len(clusters1)

    Rg = np.zeros((ny, nx))
    Rl = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):

            fname = 'R-'+str('%.6f' % G[i])+'-'+str('%.6f' % mu[j])+'.txt'
            if num_sim > 1:
                c = np.loadtxt('text/'+fname)
                c = np.mean(c, axis=0)

            else:
                c = np.loadtxt('text/'+fname)
            Rg[i, j] = c[2]
            rl = c[3:]
            Rl[i, j] = np.sum(prob*rl)

    np.savez('npz/R', Rg=Rg, Rl=Rl, g=G, mu=mu)

    if len(G) > 1:
        plot_phase_space(Rg, mu, G, 'Rg', xlabel='mu', ylabel='g',
                         xtickstep=5, ytickstep=10)
        plot_phase_space(Rl, mu, G, 'Rl', xlabel='mu', ylabel='g',
                         xtickstep=5, ytickstep=10)
    else:
        fig, ax = pl.subplots(1, figsize=(6, 4))
        ax.plot(mu, Rg[0, :], lw=2, c='r', label=r'$R_g$')
        ax.plot(mu, Rl[0, :], lw=2, c='b', label='$R_l$')
        ax.set_xlabel(r"$\omega_0$")
        ax.set_ylabel('R')

        ax.legend()
        pl.savefig('fig/'+str('%.2f.png' % G[0]), dpi=300)


# ---------------------------------------------------------------------- #
if __name__ == "__main__":

    start = time()

    adj = np.loadtxt("networks/C.txt")
    os.chdir('../data/')

    plot_correlations(G, mu)
    # plot_r_local_global_modular(G, mu, clusters1)
    plot_r_local_global_hmn(G, mu, clusters1, clusters2)
    plot_nmi(adj, G, mu, to_npz=True)

    lib.display_time(time()-start)
