#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import click
import pandas as pd
import ursgal
from scipy.stats import norm, normaltest

def normal_check(data, cutoff=0.01):
    _, p = normaltest(data)
    if p < cutoff:
        return False, p
    else:
        return True, p


def remove_outliers(data, m=8):
    median_dist = np.abs(data - np.median(data))
    dev = np.median(median_dist)
    if dev:
        s = median_dist/dev
    else:
        s = 0
    return data[s<m]


def fit_gaussian(data, file):
    data_fil = remove_outliers(data)
    is_normal, p = normal_check(data_fil)
    is_normal = True
    if is_normal:
        print(f'Normal distributed with p-value: {p}')
    else:
        print(f'Not normal distributed with p-value: {p}')
        data_fil.hist(bins=50)
        fname = os.path.splitext(os.path.basename(file))[0]
        plt.savefig(f'Offset_distribution_{fname}.png')
        exit(1)
    print(data.describe())
    print()
    print(data_fil.describe())
    mean, std = norm.fit(data_fil)
    print(f'Outlier removev: {mean}, {std}')
    return mean

@click.command()
@click.argument('files', nargs=-1)
@click.argument('database')
def main(files, database, ppm_start=-5, ppm_stop=5, no_vals=11):
    engine = 'msfragger_20190222'
    params = {
        'use_pyqms_for_mz_calculation': True,
        'database': database,
        'scan_skip_modulo_step': 2,
        'precursor_true_tolerance': 5,
        'precursor_mass_tolerance_unit': 'ppm',

        'precursor_mass_tolerance_plus': 100,
        
        'precursor_min_charge': 2,
        'precursor_max_charge': 3,
        'frag_mass_tolerance_unit': 'ppm',
        'frag_mass_tolerance': 20,
        'csv_filter_rules': [
            ['Is decoy', 'equals', 'false'],
            ['q-value', 'lte', 0.01],
        ],
        'modifications': [
            'M,opt,any,Oxidation',
            '*,opt,Prot-N-term,Acetyl',
            'C,fix,any,Carbamidomethyl',
        ]
    }
    uc = ursgal.UController(params=params, profile='QExactive+', verbose=False)
    means = []
    for file in files:
        search = uc.search(file, engine=engine, force=False)
        validated = uc.validate(search, engine='percolator_2_08')
        filtered = uc.execute_misc_engine(validated, engine='filter_csv_1_0_0')
        df = pd.read_csv(filtered)
        print('\n\n\n\n\n')
        mean = fit_gaussian(df['Accuracy (ppm)'], file)
        print(f'\n\n\n\n\n')
        means.append(mean)
    df = pd.DataFrame({'file': [os.path.basename(f) for f in files], 'offset': means})
    df.to_csv('optimal_offsets.csv', index=False)

if __name__ == '__main__':
    main()

