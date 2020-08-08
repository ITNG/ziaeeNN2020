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


def toRad(x):
    return x*(2.*np.pi)/1000.0


def toHz(x):
    return x/(2.0*np.pi)*1000


# c = np.load("R.npz")
# Rg = c['Rg'][0]
# Rl1 = c['Rl1'][0]
# Rl2 = c['Rl2'][0]
# g = c['g']
# mu = c['mu']
c = np.load("data.npz")
mu = c["mu"]
Rg = c['R_glob']
Rl1 = np.mean(c['R_loc1'], axis=0)
Rl2 = np.mean(c["R_loc2"], axis=0)
nmi1, nmi2 = c['nmi_mean']

con = np.genfromtxt("con.txt")
cor1 = np.genfromtxt("cor1.txt")
cor2 = np.genfromtxt("cor2.txt")
cor3 = np.genfromtxt("cor3.txt")
cor4 = np.genfromtxt("cor4.txt")

f = pl.figure(figsize=(8, 10))

left = 0.1
right = 0.97
fontsize = 14
labelsize = 14

gs1 = GridSpec(1, 3)
gs1.update(left=left, right=right-0.025, wspace=0.2, hspace=0.05,
           top=0.97, bottom=0.80)
ax1 = pl.subplot(gs1[0])
ax2 = pl.subplot(gs1[1])
ax3 = pl.subplot(gs1[2])

gs2 = GridSpec(2, 1)
gs2.update(left=left, right=right-0.025, hspace=0.05, top=0.77, bottom=0.3)
ax4 = pl.subplot(gs2[0])
ax5 = pl.subplot(gs2[1:])

gs3 = GridSpec(1, 4)
gs3.update(left=left, right=right-0.025, wspace=0.25, hspace=0.02,
           top=0.25, bottom=0.03)
ax6 = pl.subplot(gs3[0])
ax7 = pl.subplot(gs3[1])
ax8 = pl.subplot(gs3[2])
ax9 = pl.subplot(gs3[3])

mygraph = Image.open("graph.png")
im0 = ax1.imshow(mygraph, interpolation='none',
                 cmap='afmhot', aspect='auto')
ax1.axis('off')
ax1.text(-0.25, 0.85, "(a)", fontsize=fontsize+2, transform=ax1.transAxes)
# -------------------------------------------------------------#
im2 = ax2.imshow(con, interpolation='nearest',
                 cmap='afmhot', aspect='auto', vmax=6, vmin=0)
ax2.invert_yaxis()
cax = ax2.inset_axes([1, 0.0, 0.05, 0.4])
cbar = pl.colorbar(im2, cax=cax, ticks=[0, 3, 6])
cbar.ax.tick_params(labelsize=labelsize-2)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.text(-0.25, 0.85, "(b)", fontsize=fontsize+2, transform=ax2.transAxes)
# -------------------------------------------------------------#

mygraph = Image.open("dend.png")
im3 = ax3.imshow(mygraph, interpolation='none',
                 cmap='afmhot', aspect='auto')
ax3.axis('off')
ax3.text(-0.2, 0.85, "(c)", fontsize=fontsize+2, transform=ax3.transAxes)
ax3.text(0.35, -0.15, "Area", fontsize=fontsize, transform=ax3.transAxes)
# -------------------------------------------------------------#
ax4.plot(mu, Rg, lw=2, c='k', label=r"$R_{global}$, K=0.1")
ax4.plot(mu, Rl1, lw=2, c='r', label=r"$R_{level_1}$,  K=0.1")
ax4.plot(mu, Rl2, lw=2, c='b', label=r"$R_{level_2}$,  K=0.1")

ax4.legend(frameon=False, loc='center', fontsize=15)
# xlabels = [50, 100, 150]
# ax4.set_xticks([toRad(i) for i in xlabels])
# ax4.set_xticklabels(xlabels)
ax4.tick_params(labelsize=labelsize)
ax4.set_xticks([])
ax4.set_ylabel("R", fontsize=fontsize, rotation=0)
ax4.set_yticks([0, 0.5, 1])
ax4.set_xlim(0.1, 0.95)
ax4.set_ylim(0, 1.1)
ax4.text(-0.11, 0.85, "(d)", fontsize=fontsize+2, transform=ax4.transAxes)


idx1_ = find_first_value(mu, 0.27, tol=0.05)
idx2_ = find_first_value(mu, 0.91, tol=0.05)
ax4.fill_between(mu[idx1_:idx2_+1],
                 0, Rl1[idx1_:idx2_+1],
                 color='gray', alpha=0.2, hatch='//')

idx1 = find_first_value(mu, 0.27, tol=0.05)
idx2 = find_first_value(mu, 0.46, tol=0.05)
ax4.fill_between(mu[idx1:idx2+1],
                 0, Rl2[idx1:idx2+1],
                 color='green', alpha=0.2, hatch='\\')


# -------------------------------------------------------------#

ax5.plot(mu, nmi1, lw=2, c='k', label=r"$level_1$")
ax5.plot(mu, nmi2, lw=2, c='gray', ls='--', label=r"$level_2$")
ax5.set_ylabel('NMI', rotation=0, fontsize=fontsize+2, labelpad=5)
ax5.set_yticks([0, 0.5, 1])
xlabels = [50, 100, 150]
ax5.set_xticks([toRad(i) for i in xlabels])
ax5.set_xticklabels(xlabels)
ax5.set_xlim(0.1, 0.95)
ax5.legend(frameon=False, fontsize=15, loc='upper right')
ax5.set_xlabel(r"$\nu_0$ (Hz)", fontsize=fontsize+2, labelpad=0)
ax5.tick_params(labelsize=labelsize)
ax5.text(-0.11, 0.85, "(e)", fontsize=fontsize+2, transform=ax5.transAxes)
ax5.set_ylim(0, 1.1)


ax5.fill_between(mu[idx1_:idx2_+1],
                 0, nmi1[idx1_:idx2_+1],
                 color='gray', alpha=0.2, hatch='//')

ax5.fill_between(mu[idx1:idx2+1],
                 0, nmi2[idx1:idx2+1],
                 color='green', alpha=0.2, hatch='\\')

# -------------------------------------------------------------#


def cor_plot(ax, cor, vmax=None, vmin=None, lb=None,
             ticks=None, format='%.0f', ylabel=False):

    im = ax.imshow(cor,
                   interpolation='nearest',
                   cmap='afmhot',
                   alpha=0.9,
                   vmin=vmin, vmax=vmax,
                   )
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    divider = make_axes_locatable(ax)
    # cax = ax.inset_axes([1, 0.0, 0.05, 0.3])
    # cbar = pl.colorbar(im, cax=cax, ticks=ticks, format=format)
    cax = divider.append_axes("top", size="5%", pad=0.05)
    cbar = pl.colorbar(im, cax=cax,
                       ticks=ticks,
                       orientation='horizontal',
                       format=format)
    cbar.ax.tick_params(labelsize=labelsize-2)
    cbar.ax.xaxis.set_ticks_position('top')
    ax.text(-0.22, 0.85, lb,
            fontsize=labelsize+2,
            transform=ax.transAxes)
    ax.set_xlabel("index", fontsize=14)
    if ylabel:
        ax.set_ylabel("index", fontsize=14)


cor_plot(ax6, cor1, vmax=1, vmin=0.97, lb='(f)', ylabel=True,
         ticks=[0.97, 1], format="%.2f")
cor_plot(ax7, cor2, vmax=1, vmin=-1, lb='(g)',
         ticks=[-1, 1])
cor_plot(ax8, cor3, vmax=1, vmin=-1, lb='(h)',
         ticks=[-1, 1])
cor_plot(ax9, cor4, vmax=1, vmin=-1, lb='(i)',
         ticks=[-1, 1])

pl.savefig('../000.pdf', dpi=150)
