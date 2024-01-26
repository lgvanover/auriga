'''
loads the evolution of median epsilon for one specific halo

in
'''
import sys
sys.path.insert(0, '/home/lgvanover/intro_scripts/')
import pandas as pd
import matplotlib.pyplot as plt

def create_plot(halo, timescale, level2=False, level3=False, level4=False):

    levels = []
    if level2:
        levels.append(2)
    if level3:
        levels.append(3)
    if level4:
        levels.append(4)

    for n in levels:
        df = pd.read_hdf("/home/lgvanover/intro_scripts/halos/new/level"+str(x)+"_halo")

    df_level3 = pd.read_hdf("halos/new/level3_halo6_all_new.hdf5", "table")
    df_level4 = pd.read_hdf("halos/level4/halo_6.hdf5", "table")
    df_level2 = pd.read_hdf("halos/new/level2_halo6_all_new.hdf5", "table")

    plt.plot(df_level2['Redshift'], df_level2['Epsilon'], 'r', label='Level 2')
    plt.plot(df_level3['Redshift'], df_level3['Epsilon'], 'b', label='Level 3')
    plt.plot(df_level4['Redshift'], df_level4['Med_Eps'], 'g', label='Level 4')

    xlim(6,0)
    ylim(-.1,1.0)
    plt.legend(loc=0)
    plt.xlabel("Redshift")
    plt.ylabel("Median Epsilon")
    plt.title("Halo 6")