import sys
from loadmodules import *
import pandas as pd
# matplotlib.use("pdf")
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     # 'font.family': 'serif',
#     # 'text.usetex': True,
#     # 'pgf.rcfonts': False,
# })


target_gas_mass = 3.74806e-06
hubbleparam = 0.67769999999999997
snap = 127
nhalos = 30

import shutil
sys.path.insert(0, '../..')
# src = '/cds/kravtsov/Simulations/Auriga/satellite_utilities.py'
# dst =  '/home/lgvanover/intro_scripts/satellite_utilities.py'
# shutil.copyfile(src, dst)
import satellite_utilities as sat


dpp = loadDefaultPlotParams('mnras')
print(dpp)
figwidth = dpp['pagewidth']
print(figwidth)

rcParams["text.usetex"] = False

# plt.style.use(dpp)

fig, axs = plt.subplots(2,4,figsize=(figwidth,.7*figwidth))
axs = axs.ravel()

colors = plt.cm.plasma(np.linspace(0, 1, 41))
for i, snap in enumerate(range(22,62)):

    df = pd.read_csv("/home/lgvanover/intro_scripts/halos/new/distr/halo6_snap"+str(snap)+".csv")
    edges = df['Bin_Edges'][:-1]
    y = df['Bin_Y'][:-1]
    axs[1].plot(edges, y, color=colors[i])

    dfyoung = pd.read_hdf("/home/lgvanover/intro_scripts/halos/new/distr/age_filt/age10/halo6_snap"+str(snap)+".hdf5", 'table')
    edgesyoung = dfyoung['Bin_Edges'][:-1]
    yyoung = dfyoung['Bin_Y'][:-1]
    axs[2].plot(edgesyoung, yyoung, color=colors[i])

    dfyoungest = pd.read_hdf("/home/lgvanover/intro_scripts/halos/new/distr/age_filt/age1/halo6_snap"+str(snap)+".hdf5", 'table')
    edgesyoungest = dfyoungest['Bin_Edges'][:-1]
    yyoungest = dfyoungest['Bin_Y'][:-1]
    axs[3].plot(edgesyoungest, yyoungest, color=colors[i])

    # now for halo 24 
    df = pd.read_csv("/home/lgvanover/intro_scripts/halos/new/distr/halo24_snap"+str(snap)+".csv")
    edges = df['Bin_Edges'][:-1]
    y = df['Bin_Y'][:-1]
    axs[5].plot(edges, y, color=colors[i])

    dfyoung = pd.read_hdf("/home/lgvanover/intro_scripts/halos/new/distr/age_filt/age10/halo24_snap"+str(snap)+".hdf5",'table')
    edgesyoung = dfyoung['Bin_Edges'][:-1]
    yyoung = dfyoung['Bin_Y'][:-1]
    axs[6].plot(edgesyoung, yyoung, color=colors[i])

    dfyoungest = pd.read_hdf("/home/lgvanover/intro_scripts/halos/new/distr/age_filt/age1/halo24_snap"+str(snap)+".hdf5", 'table')
    edgesyoungest = dfyoungest['Bin_Edges'][:-1]
    yyoungest = dfyoungest['Bin_Y'][:-1]
    axs[7].plot(edgesyoungest, yyoungest, color=colors[i])

df_6 = pd.read_hdf("/home/lgvanover/intro_scripts/halos/new/level3_halo6_all_new.hdf5", "table")
z6 = df_6['Redshift']
lookbacktimes6 = sat.return_lookbacktime_from_a((z6+1.0)**(-1.0))
axs[0].plot(lookbacktimes6,df_6['Epsilon'],label='Median Eps')

df_24 = pd.read_hdf("/home/lgvanover/intro_scripts/halos/new/level3_halo24_all_new.hdf5", "table")
z24 = df_24['Redshift']
lookbacktimes24 = sat.return_lookbacktime_from_a((z24+1.0)**(-1.0))
axs[4].plot(lookbacktimes24,df_24['Epsilon'],label='Median Eps')

for x in [0, 4]:
    axs[x].set_xlim(12.5, 0)
    axs[x].set_ylim(-.2, 1)
    if x < 4:
        axs[x].get_xaxis().set_ticklabels([])

for x in [1, 2, 3, 5, 6, 7]:
    axs[x].set_xlim(-2, 2)
    axs[x].set_ylim(0, 3.5)
    if x < 4:
        axs[x].get_xaxis().set_ticklabels([])

axs[0].set_ylabel('Median $\epsilon$')
axs[4].set_ylabel('Fraction of Starts')
axs[1].set_title('No Age Filter', x=0.65, y=.85, fontsize=8)
axs[2].set_title('$<$ 1 Gyr', x=0.75, y=.85, fontsize=8)
axs[3].set_title('$<$ .1 Gyr', x=0.75, y=.85, fontsize=8)
# axs[3].set_title('$<$ .1 Gyr', x=0.75, y=.85)


axs[0].set_xlabel('Halo 6', x=.75, y=.2)
axs[0].xaxis.set_label_coords(0.75,0.2)
axs[4].set_title('Halo 24', x=.75, y=.2, fontsize=8)
# axs[4].set_xlabel('Halo 24', fontsize=7, rotation=0)
# axs[4].yaxis.set_label_coords(0.75,0.2)

axs[2].yaxis.set_tick_params(labelleft=False)
axs[3].yaxis.set_tick_params(labelleft=False)
axs[6].yaxis.set_tick_params(labelleft=False)
axs[7].yaxis.set_tick_params(labelleft=False)

for index, label in enumerate(axs[5].xaxis.get_ticklabels()):
    if index % 2 != 0:
        label.set_visible(False)

for index, label in enumerate(axs[6].xaxis.get_ticklabels()):
    if index % 2 != 0:
        label.set_visible(False)

for index, label in enumerate(axs[7].xaxis.get_ticklabels()):
    # print(index, label)
    if index % 2 != 0:
        label.set_visible(False)

axs[4].set_xlabel('Look Back Time (Gyrs)')
axs[5].set_xlabel("$\epsilon$")
axs[6].set_xlabel("$\epsilon$")
axs[7].set_xlabel("$\epsilon$")

xlim(-2,2)
# fig.supylabel("Number")
# fig.supxlabel("$\epsilon$")
# fig.suptitle("Distribution of Circularity")
# fig.subplots_adjust(hspace = 0.1, wspace=0.1)

# print("hello")
# # plt.show()
# print(rcParams["text.usetex"])

fig.tight_layout()
fig.savefig("/home/lgvanover/intro_scripts/figures/figs/distr_eps.pdf", format="pdf")
plt.show()