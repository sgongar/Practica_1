# import millennium_query
from os import getcwd
import pandas as pd

from math import log
from os import getcwd

import numpy as np

import matplotlib.pyplot as plt


def output_galaxies(galaxy_classification, label, sizes):
    """

    """
    data_dir = f'{getcwd()}/datos/{galaxy_classification}.csv'
    data = pd.read_csv(data_dir)

    fig, ax = plt.subplots(figsize=(11.69, 8.27))

    redshift_list = []
    mass_list = []
    final_mass = 0.0

    try:
        for idx, row in data.iterrows():
          # Check if we are at first position
            if int(row['snapnum']) == 63:
                if len(redshift_list) != 0:
                    ax.plot(redshift_list, mass_list)
            redshift_list = []
            mass_list = []
            if float(row['stellarMass']) != 0.0:
                final_mass = float(row['stellarMass'])
            else:
                raise Exception
            if final_mass != 0.0:
                redshift_list.append(log(1+row['redshift']))
                mass_list.append(log(row['stellarMass']/final_mass))
        # We are not at the first element
        else:
            if row['stellarMass'] != 0.0:
                redshift_list.append(log(1+row['redshift']))
                mass_list.append(log(row['stellarMass']/final_mass))
    except Exception as e:
        print(str(e))

    ax.set_title(label)

    ax.set_xlabel(r'$\log [1+Z]$')
    ax.set_ylabel(r'$\log[M_{main} (z) / M_0]$')

    ax.grid(True)
    plt.savefig(f'{galaxy_classification}.png', dpi=72, quality=80,
                optimize=True, progressive=True)  


def output_halos(halo_classification, label, sizes):
    """

    """

    data_dir = f'{getcwd()}/datos/{halo_classification}.csv'
    data = pd.read_csv(data_dir)

    fig, ax = plt.subplots(figsize=(11.69, 8.27))

    redshift_list = []
    mass_list = []
    last_progenitor_id_list = []
    final_mass = 0.0
    for idx, row in data.iterrows():
        # Check if we are at first position
        if int(row['snapNum']) == 63:
            if len(redshift_list) != 0:
                print(last_progenitor_id)
                ax.plot(redshift_list, mass_list)
            redshift_list = []
            mass_list = []
            if float(row['m_Mean200']) != 0.0:
                final_mass = float(row['m_Mean200'])
            else:
                raise Exception
            if final_mass != 0.0:
                redshift_list.append(log(1+row['redshift']))
                mass_list.append(log(row['m_Mean200']/final_mass))
                last_progenitor_id = row['lastProgenitorId']
        # We are not at the first element
        else:
            if row['m_Mean200'] != 0.0:
                redshift_list.append(log(1+row['redshift']))
                mass_list.append(log(row['m_Mean200']/final_mass))

    ax.set_title(label)
    y_ticks = np.arange(sizes['y'][0], sizes['y'][1], sizes['y'][2])
    x_ticks = np.arange(sizes['x'][0], sizes['x'][1], sizes['x'][2])

    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xlabel(r'$\log [1+Z]$')
    ax.set_ylabel(r'$\log[M_{main} (z) / M_0]$')

    ax.grid(True)
    plt.savefig(f'{halo_classification}.png', dpi=72, quality=80,
                optimize=True, progressive=True)


# Rutina principal para crear las gráficas
def main():
    # Nombre de los ficheros de los halos
    halos_classification = ['halos_10', 'halos_10_11','halos_11_12',
                            'halos_12_13', 'halos_13']

    # Etiquetas de las graficas de los halos
    halos_label = [r'Halos < ' + r'$10^{10}$',
                r'$10^{10}$' + ' < Halos < ' + r'$10^{11}$',
                r'$10^{11}$' + ' < Halos < ' + r'$10^{12}$',
                r'$10^{12}$' + ' < Halos < ' + r'$10^{13}$',
                r'Halos > ' + r'$10^{13}$',]

    # Tamaño de las graficas de los halos
    halos_size = [{'x': [0, 0.75, 0.25], 'y': [-1, 1, 0.5]},
                {'x': [0, 2.0, 0.25], 'y': [-2.5, 1, 0.5]},
                {'x': [0, 2.75, 0.25], 'y': [-6.5, 1, 0.5]},
                {'x': [0, 2.75, 0.25], 'y': [-7.5, 1, 0.5]},
                {'x': [0, 2.75, 0.25], 'y': [-10, 1, 0.5]}]

    # Bucle for para crear las graficas de los halos
    for idx, halo_classification in enumerate(halos_classification):
        output_halos(halo_classification, halos_label[idx], halos_size[idx])

    # Nombre de los ficheros de las galaxias
    galaxies_classification = ['galaxias_10_11', 'galaxias_11_12',
                               'galaxias_12_13', 'galaxias_13']

    # Etiquetas de las gráficas de las galaxias
    galaxies_label = [r'$10^{10}$' + ' < Galaxias < ' + r'$10^{11}$',
                      r'$10^{11}$' + ' < Galaxias < ' + r'$10^{12}$',
                      r'$10^{12}$' + ' < Galaxias < ' + r'$10^{13}$',
                      r'Galaxias > ' + r'$10^{13}$']

    # Tamaño de las gráficas de las galaxias
    galaxies_size = [{'x': [0, 2.75, 0.25], 'y': [-9.5, 1, 0.5]},
                    {'x': [0, 2.75, 0.25], 'y': [-9.5, 1, 0.5]},
                    {'x': [0, 2.75, 0.25], 'y': [-9.5, 1, 0.5]},
                    {'x': [0, 2.75, 0.25], 'y': [-12, 1, 0.5]}]

    # Bucle for para crear las graficas de las galaxias
    try:
        for idx, galaxy_classification in enumerate(galaxies_classification):
            output_galaxies(galaxy_classification, galaxies_label[idx], 
            galaxies_size[idx])
    except Exception as e:
        print(str(e))

if __name__ == "__main__":
    main()
