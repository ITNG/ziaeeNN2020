import os
import lib
import matplotlib
import numpy as np
import pylab as pl
from PIL import Image
from PIL import ImageOps
from cycler import cycler
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

os.chdir('data/')

# ------------------------------------------------------- #

f = pl.figure(figsize=(8, 10))

left = 0.1
right = 0.97
fontsize = 14
labelsize = 14


gs1 = GridSpec(1, 3)
gs1.update(left=left, right=right, wspace=0.1, hspace=0.05,
           top=0.97, bottom=0.80)
ax1 = pl.subplot(gs1[0])
ax2 = pl.subplot(gs1[1])
ax3 = pl.subplot(gs1[2])

gs2 = GridSpec(2, 1)
gs2.update(left=left, right=right-0.028, hspace=0.05, top=0.77, bottom=0.42)
ax4 = pl.subplot(gs2[0])
ax5 = pl.subplot(gs2[1:])

gs3 = GridSpec(1, 3)
gs3.update(left=left, right=right-0.028, wspace=0.4, hspace=0.02,
           top=0.37, bottom=0.2)
ax6 = pl.subplot(gs3[0])
ax7 = pl.subplot(gs3[1])
ax8 = pl.subplot(gs3[2])

gs4 = GridSpec(1, 3)
gs4.update(left=left, right=right-0.028, wspace=0.4,
           top=0.175, bottom=0.05)
ax9 = pl.subplot(gs4[0])
ax10 = pl.subplot(gs4[1])
ax11 = pl.subplot(gs4[2])


communities = lib.read_from_file('communities.txt', 5)
con = np.loadtxt("r_C65.dat")
length = np.loadtxt("r_L65.dat")

c = np.load("R.npz")
Rg = c['Rg'][0]
Rl1 = c['Rl1'][0]
Rl2 = c['Rl2'][0]
g = c['g']
mu = c['mu']

cor1 = np.load("cor-0.06.npz")['cor']
cor2 = np.load("cor-0.29.npz")['cor']
cor3 = np.load("cor-0.63.npz")['cor']

c1 = np.load('nmi-30.npz')
nmi = c1['nmi']
MU = c1['omega']
std = c1["std"]


mygraph = Image.open("graph.png")
border = (1, 1, 1, 1)  # left, up, right, bottom
cropped_graph = ImageOps.crop(mygraph, border)
im0 = ax1.imshow(cropped_graph, interpolation='none',
                 cmap='afmhot', aspect='equal')
ax1.axis('off')
ax1.text(-0.25, 0.85, "(a)", fontsize=fontsize+2, transform=ax1.transAxes)
# -------------------------------------------------------------#
lib.plot_rectangles(communities, np.log10(con), ax2, ylabel=True,
                    title="Coupling weight")

ax2.text(-0.25, 0.85, "(b)", fontsize=fontsize+2, transform=ax2.transAxes)
# -------------------------------------------------------------#
lib.plot_rectangles(communities, np.log10(length), ax3,
                    title="Distance")

ax3.text(-0.2, 0.85, "(c)", fontsize=fontsize+2, transform=ax3.transAxes)
# -------------------------------------------------------------#
ax4.plot(mu, Rg, lw=2, c='k', label=r"$R_{global}$")
ax4.plot(mu, Rl1, lw=2, c='r', label=r"$R_{level_1}$")
ax4.plot(mu, Rl2, lw=2, c='b', label=r"$R_{level_2}$")

ax4.legend(frameon=False, loc='upper right', fontsize=15)
xlabels = [50, 100, 150]
ax4.set_xticks([lib.toRad(i) for i in xlabels])
ax4.set_xticklabels(xlabels)
ax4.tick_params(labelsize=labelsize)
ax4.set_xticks([])
ax4.set_ylabel("R", fontsize=fontsize, rotation=0)
ax4.set_yticks([0, 0.5, 1])
# ax4.set_xlim(0.04, 0.85)
ax4.set_xlim(0.0, 0.99)
ax4.set_ylim(0, 1.1)
ax4.text(-0.11, 0.85, "(d)", fontsize=fontsize+2, transform=ax4.transAxes)
# -------------------------------------------------------------#

ax5.errorbar(MU, nmi[0], yerr=std[0], c="k",
             elinewidth=1,
             capsize=2,
             label=r"$level_1$")
ax5.errorbar(MU, nmi[1], yerr=std[1], c="gray",
             elinewidth=1,
             capsize=2,
             label=r"$level_2$")

ax5.set_ylabel('NMI', rotation=0, fontsize=fontsize+2, labelpad=5)
ax5.set_yticks([0, 0.2, 0.4])
xlabels = [0, 25, 50, 75, 100, 125, 150]
ax5.set_xticks([lib.toRad(i) for i in xlabels])
ax5.set_xticklabels(xlabels)
# ax5.set_xlim(0.04, 0.85)
ax5.set_xlim(0.0, 0.99)
ax5.legend(frameon=False, fontsize=15, loc='upper right')
ax5.set_xlabel(r"$\nu_0$ (Hz)", fontsize=fontsize+2, labelpad=0)
ax5.tick_params(labelsize=labelsize)
ax5.text(-0.11, 0.85, "(e)", fontsize=fontsize+2, transform=ax5.transAxes)
ax5.set_ylim(0, 0.5)
# -------------------------------------------------------------#

lib.cor_plot(ax6, cor1, vmax=1, vmin=-1, lb='(f)', ticks=[-1, 1], ylabel=True)
lib.cor_plot(ax7, cor2, vmax=1, vmin=-1, lb='(g)', ticks=[-1, 1])
lib.cor_plot(ax8, cor3, vmax=1, vmin=-1, lb='(h)', ticks=[-1, 1])

logscale = True
clusters1 = [12, 13, 14, 14, 12]
clusters2 = [32, 33]

ylim = (5e-4, 3e-1)
yticks = [1e-3, 1e-2, 1e-1]

lib.plot_correlation_distribution(
    cor1,
    clusters1,
    ax9,
    nbins=50,
    xticks=[0, 1],
    yticks=yticks,
    ylim=ylim,
    bin_range=(0, 1),
    xlim=(0, 1.01),
    label=0.06,
    xlabel=True,
    ylabel=True,
    logscale=logscale,
    label1="(i)"
)
lib.plot_correlation_distribution(
    cor2,
    clusters1,
    ax10,
    yticks=yticks,
    ylim=ylim,
    nbins=50,
    xticks=[-1, 0, 1],
    xlim=(-1, 1),
    # bin_range=(-0.5, 1),
    label=0.29,
    xlabel=True,
    logscale=logscale,
    label1="(j)"
)


lib.plot_correlation_distribution(
    cor3,
    clusters1,
    ax11,
    nbins=50,
    yticks=yticks,
    ylim=ylim,
    # bin_range=(-0.5, 0.5),
    xlim=(-1, 1),
    xticks=[-1, 0, 1],
    label=0.63,
    xlabel=True,
    logscale=logscale,
    label1="(k)"
)
pl.tight_layout()
pl.savefig('../004.pdf', dpi=150)
