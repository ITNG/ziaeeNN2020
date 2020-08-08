import os
import lib
import matplotlib
import numpy as np
import pylab as pl
from PIL import Image
from cycler import cycler
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

os.chdir('data/')


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx  # , array[idx]


def find_first_value(arr, value, tol=0.01):
    for i, item in enumerate(arr):
        if np.abs(item-value) < tol:
            return i
    print ("match not found!")
    exit(0)


def fill_between(Rl, ax, thr1=0.95, thr2=0.6, tol=0.02,
                 color='gray', pattern='//'):

    indecies = []
    idx1 = find_first_value(Rl, thr1, tol=tol)
    idx2 = find_first_value(Rl, thr2, tol=tol)
    ax.fill_between(mu[idx1:idx2+1], 0, Rl[idx1:idx2+1],
                    color=color, alpha=0.2, hatch=pattern)
    indecies.append(idx1)
    indecies.append(idx2+1)

    return indecies


def make_ticklabels_invisible(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_visible(False)


def toRad(x):
    return x*(2.*np.pi)/1000.0


def toHz(x):
    return x/(2.0*np.pi)*1000


c1 = np.load("data1.npz")
Rg1 = c1['R_glob']
Rl1 = np.mean(c1['R_loc'], axis=0)
nmi1 = c1['nmi_mean']
c2 = np.load("data2.npz")
Rg2 = c2['R_glob']
Rl2 = np.mean(c2['R_loc'], axis=0)
nmi2 = c2["nmi_mean"]
g = c1['g']
mu = c1['mu']

con = np.genfromtxt("C.dat")
cor1 = np.load("cor1.npz")['cor']
cor2 = np.load("cor2.npz")['cor']
cor3 = np.load("cor3.npz")['cor']

# nmi = np.loadtxt('nmi.txt')

d2 = [2, 3.2, 4.4, 5.6]

f = pl.figure(figsize=(8, 12))

left = 0.1
right = 0.97
fontsize = 14
labelsize = 14

gs1 = GridSpec(1, 3)
gs1.update(left=left, right=right-0.025, wspace=0.4, hspace=0.05,
           top=0.97, bottom=0.85)
ax1 = pl.subplot(gs1[0])
ax2 = pl.subplot(gs1[1])
ax3 = pl.subplot(gs1[2])

gs2 = GridSpec(3, 1)
gs2.update(left=left, right=right-0.025, hspace=0.05, top=0.81, bottom=0.4)
ax4 = pl.subplot(gs2[0])
ax5 = pl.subplot(gs2[1])
ax6 = pl.subplot(gs2[2])

gs3 = GridSpec(1, 3)
gs3.update(left=left, right=right-0.025, wspace=0.32, hspace=0.02,
           top=0.35, bottom=0.2)
ax7 = pl.subplot(gs3[0])
ax8 = pl.subplot(gs3[1])
ax9 = pl.subplot(gs3[2])

gs4 = GridSpec(1, 3)
gs4.update(left=left, right=right-0.025, wspace=0.3,
           top=0.17, bottom=0.05)
ax10 = pl.subplot(gs4[0])
ax11 = pl.subplot(gs4[1])
ax12 = pl.subplot(gs4[2])


mygraph = Image.open("graph.png")
im0 = ax1.imshow(mygraph, interpolation='none',
                 cmap='afmhot', aspect='auto')  # , cmap=pl.cm.ocean
ax1.axis('off')
ax1.text(-0.25, 0.85, "(a)", fontsize=fontsize+2, transform=ax1.transAxes)
# -------------------------------------------------------------#
im2 = ax2.imshow(con, interpolation='nearest',
                 cmap='afmhot', aspect='auto')  # , cmap=pl.cm.ocean
ax2.invert_yaxis()
cax = ax2.inset_axes([1, 0.0, 0.05, 0.3])
cbar = pl.colorbar(im2, cax=cax, ticks=[0, 1])
cbar.ax.tick_params(labelsize=labelsize)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.text(-0.25, 0.85, "(b)", fontsize=fontsize+2, transform=ax2.transAxes)
# -------------------------------------------------------------#
NUM_COLORS = 4
cm = pl.get_cmap('gist_rainbow')  # , gist_yarg
ax3.set_prop_cycle(cycler('color', [cm(1.*ii/NUM_COLORS)
                                    for ii in range(NUM_COLORS)]))
labels = ["2.0", "3.2", "4.4", "5.6"]
for i in range(4):
    nmi = np.load("nmi/"+labels[i]+".npz")
    ax3.plot(nmi['mu'], nmi["nmi"], lw=2, label=str(d2[i]))

ax3.set_xlabel(r"$\nu_0$ (Hz)", fontsize=fontsize+1, labelpad=-10)
ax3.set_ylabel("NMI", fontsize=fontsize-1, rotation=0, labelpad=-10)
ax3.set_yticks([0.0, 0.4, 0.8])
ax3.set_ylim(0, 1.1)
ax3.legend(frameon=False, loc='upper right', fontsize=14)
ax3.set_xlim(0.11, 0.99)
ax3.set_ylim(0, 1)
xlabels = [50, 150]
ax3.set_xticks([toRad(i) for i in xlabels])
ax3.set_xticklabels(xlabels)
ax3.text(-0.2, 0.92, "(c)", fontsize=fontsize+2, transform=ax3.transAxes)
ax3.tick_params(labelsize=labelsize)
# -------------------------------------------------------------#

ax4.plot(mu, Rg1, lw=2, c='k', label=r"$R_{global}$, K=0.1")
ax4.plot(mu, Rl1, lw=2, c='r', label=r"$R_{local}$,   K=0.1")
ax4.legend(frameon=False, loc='upper right', fontsize=15)
xlabels = [50, 100, 150]
ax4.set_xticks([toRad(i) for i in xlabels])
ax4.set_xticklabels(xlabels)
ax4.tick_params(labelsize=labelsize)
ax4.set_ylabel("R", fontsize=fontsize, rotation=0)
ax4.set_yticks([0, 0.5, 1])
ax4.set_xlim(0.0, 0.95)
ax4.set_ylim(0, 1.1)
ax4.set_xticks([])
ax4.text(-0.11, 0.85, "(d)", fontsize=fontsize+2, transform=ax4.transAxes)


# i_nmi1 = fill_between(Rl1, ax4, tol=0.05,
#                       thr1=0.98, thr2=0.3)
idx1 = find_first_value(mu, 0.3, tol=0.05)
idx2 = find_first_value(mu, 0.835, tol=0.05)
ax4.fill_between(mu[idx1:idx2+1],
                 0, Rl1[idx1:idx2+1],
                 color='gray', alpha=0.2, hatch='//')

# -------------------------------------------------------------#

ax5.plot(mu, Rg2, lw=2, c='k', label=r"$R_{global}$, K=0.5")
ax5.plot(mu, Rl2, lw=2, c='r', label=r"$R_{local}$,   K=0.5")
ax5.legend(frameon=False, loc='upper right', fontsize=15)
xlabels = [50, 100, 150]
ax5.set_xticks([toRad(i) for i in xlabels])
ax5.set_xticklabels(xlabels)
ax5.tick_params(labelsize=labelsize)
ax5.set_ylabel("R", fontsize=fontsize, rotation=0)
ax5.set_yticks([0, 0.5, 1])
ax5.set_xlim(0.0, 0.95)
ax5.set_ylim(0, 1.1)
ax5.set_xticks([])
ax5.text(-0.11, 0.85, "(e)", fontsize=fontsize+2, transform=ax5.transAxes)

# i_nmi2 = fill_between(Rl[1, :], ax5, tol=0.02,
#                       thr1=0.965, thr2=0.7,
#                       color='green', pattern='\\')
idx1_ = find_first_value(mu, 0.385, tol=0.05)
idx2_ = find_first_value(mu, 0.885, tol=0.05)
ax5.fill_between(mu[idx1_:idx2_+1],
                 0, Rl2[idx1_:idx2_+1],
                 color='green', alpha=0.2, hatch='\\')


# -------------------------------------------------------------#

ax6.plot(mu, nmi1, lw=2, c='k', label="K=0.1")
ax6.plot(mu, nmi2, lw=2, c='gray', ls='--', label="K=0.5")
ax6.set_ylabel('NMI', rotation=0, fontsize=fontsize+2, labelpad=5)
ax6.set_yticks([0, 0.5, 1])
xlabels = [50, 100, 150]
ax6.set_xticks([toRad(i) for i in xlabels])
ax6.set_xticklabels(xlabels)
ax6.set_xlim(0.0, 0.95)
ax6.legend(frameon=False, loc='upper right', fontsize=15)
ax6.set_xlabel(r"$\nu_0$ (Hz)", fontsize=fontsize+1, labelpad=1)
ax6.tick_params(labelsize=labelsize)
ax6.text(-0.11, 0.85, "(f)", fontsize=fontsize+2, transform=ax6.transAxes)
ax6.set_ylim(0, 1.1)

ax6.fill_between(mu[idx1:idx2+1],
                 0, nmi1[idx1:idx2+1],
                 color='gray', alpha=0.2, hatch='//')

ax6.fill_between(mu[idx1_:idx2_+1],
                 0, nmi2[idx1_:idx2_+1],
                 color='green', alpha=0.2, hatch='\\')
# -------------------------------------------------------------#


def cor_plot(ax, cor, vmax=None, vmin=None, lb=None,
             ticks=None, format='%.0f', ylabel=False):

    im = ax.imshow(cor,
                   aspect='auto',
                   interpolation='nearest',
                   cmap='afmhot',
                   alpha=0.9,
                   vmin=vmin, vmax=vmax,
                   )
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    divider = make_axes_locatable(ax)
    cax = ax.inset_axes([1, 0.0, 0.05, 0.3])
    cbar = pl.colorbar(im, cax=cax, ticks=ticks)
    # cax = divider.append_axes("bottom", size="5%", pad=0.05)
    # cbar = pl.colorbar(im, cax=cax,
    #                    ticks=ticks,
    #                    orientation='horizontal',
    #                    format=format)
    cbar.ax.tick_params(labelsize=13)
    ax.text(-0.25, 0.85, lb,
            fontsize=15,
            transform=ax.transAxes)
    ax.set_xlabel("index", fontsize=14)
    if ylabel:
        ax.set_ylabel("index", fontsize=14)


cor_plot(ax7, cor1, vmax=1, vmin=0.98, lb='(g)', ylabel=True,
         ticks=[0.98, 1], format="%.2f")
cor_plot(ax8, cor2, vmax=1, vmin=-1, lb='(h)', ticks=[-1, 1])
cor_plot(ax9, cor3, vmax=1, vmin=-1, lb='(i)', ticks=[-1, 1])

clusters1 = [20]*3
logscale = True
ylim = (5e-4, 3e-1)
yticks = [1e-3, 1e-2, 1e-1]

lib.plot_correlation_distribution(
    cor1,
    clusters1,
    ax10,
    nbins=50,
    xticks=[0, 1],
    yticks=yticks,
    ylim=ylim,
    bin_range=(0., 1),
    xlim=(0, 1.01),
    label=10,
    label1="(j)",
    xlabel=True,
    ylabel=True,
    logscale=logscale,
)
lib.plot_correlation_distribution(
    cor2,
    clusters1,
    ax11,
    nbins=50,
    label=72,
    yticks=yticks,
    ylim=ylim,
    xlim=(-1, 1.01),
    label1="(k)",
    xlabel=True,
    logscale=logscale,
)
lib.plot_correlation_distribution(
    cor3,
    clusters1,
    ax12,
    nbins=50,
    xlim=(-0.5, 0.5),
    label=150,
    label1="(l)",
    xlabel=True,
    yticks=yticks,
    ylim=ylim,
    logscale=logscale,
)

pl.savefig('../002.pdf', dpi=150)
