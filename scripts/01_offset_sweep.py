#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import click
import pandas as pd
import ursgal
from scipy.stats import norm, normaltest
import sys


def draw_heatmap(plot_name, tol_dict, prec_tolerances, frag_tolerances):
    z = np.zeros(
        (len(prec_tolerances), len(frag_tolerances))
    )
    for i, prec_tol in enumerate(prec_tolerances):
        for j, frag_tol in enumerate(frag_tolerances):
            z[i, j] = tol_dict[(prec_tol, frag_tol)]
    fig, ax = plt.subplots()
    im = ax.imshow(z)
    ax.set_xticklabels(['-'] + frag_tolerances)
    ax.set_yticklabels(['-'] + prec_tolerances)
    fig.colorbar(im)
    plt.xlabel('frag tolerance [ppm]')
    plt.ylabel('prec tolerance [ppm]')
    plt.title(plot_name)
    plt.savefig(plot_name)
    return z


@click.command()
@click.argument('files', nargs=-1)
@click.argument('database')
@click.option('--force', '-f', is_flag=True, help="Print more output.")
def main(files, database, force):
    # engine = 'msgfplus_v2019_04_18'
    # engine = 'msfragger_20190222'
    engine = 'xtandem_alanine'
    params = {
        # 'use_pyqms_for_mz_calculation': True,
        'database': database,
        'scan_skip_modulo_step': 2,
        'precursor_true_tolerance': 5,
        'precursor_mass_tolerance_unit': 'ppm',
        'precursor_mass_tolerance_plus': 300,
        'precursor_mass_tolerance_minus': 300,
        # 'precursor_mass_tolerance_minus': 5,
        # 'precursor_mass_tolerance_plus': 5,

        'frag_mass_tolerance_unit': 'ppm',
        'frag_mass_tolerance': 700,
        # 'frag_mass_tolerance': 20,
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

    # frag tol 700ppm
    # prec tol 300ppm
    uc = ursgal.UController(params=params, profile='QExactive+', verbose=False)
    means = []

    # in uncomment this, remember to change frag_tolerance unit!
    # frag_tolerances_dalton = [50, 25, 17.5, 10, 5, 2.5]
    precursor_tolerances = [20, 17.5, 15  , 12.5, 10  , 5 ]
    frag_tolerances      = [25, 20  , 17.5, 15  , 12.5, 10]
    all_matrices = []
    for file in files:
        if '01' in file:
            continue
        uc.params['machine_offset_in_ppm'] = 0
        uc.params['prefix'] = ''
        uc.params['precursor_mass_tolerance_unit'] = 'ppm'
        uc.params['precursor_mass_tolerance_plus'] = 300
        uc.params['precursor_mass_tolerance_minus'] = 300
        uc.params['frag_mass_tolerance'] = 700
        uc.params['frag_mass_tolerance_unit'] = 'ppm'
        # 'precursor_mass_tolerance_minus': 5,
        # 'precursor_mass_tolerance_plus': 5,

        search = uc.search(file, engine=engine, force=force)

        validated = uc.validate(search, engine='percolator_3_4', force=force)
        # validated = uc.validate(search, engine='percolator_2_08', force=force)

        filtered = uc.execute_misc_engine(validated, engine='filter_csv_1_0_0', force=force)
        df = pd.read_csv(filtered)
        machine_offset = -1 * df['Accuracy (ppm)'].median()
        means.append(machine_offset)
        tolerance_dict = {}
        uc.params['prefix'] = 'corrected_mgf'
        uc.params['machine_offset_in_ppm'] = machine_offset
        corrected_mgf = uc.convert_to_mgf_and_update_rt_lookup(file, force=force)
        for prec_tol in precursor_tolerances:
            uc.params['precursor_mass_tolerance_plus'] = prec_tol
            uc.params['precursor_mass_tolerance_minus'] = prec_tol
            for frag_tol in frag_tolerances:
                # if prec tol > frag_tol:
                #     continue
                uc.params['prefix'] = f'corrected_mgf_{prec_tol}_{frag_tol}'
                uc.params['frag_mass_tolerance'] = frag_tol
                search = uc.search(corrected_mgf, engine=engine, force=force)
                validated = uc.validate(search, engine='percolator_3_4', force=force)
                # validated = uc.validate(search, engine='percolator_2_08', force=force)
                filtered = uc.execute_misc_engine(validated, engine='filter_csv_1_0_0', force=force)
                results = pd.read_csv(filtered)
                tolerance_dict[(prec_tol, frag_tol)] = len(results)
                file_basename = os.path.splitext(os.path.basename(file))[0]
                plot_name = f'{file_basename}_heatmap.png'
        z = draw_heatmap(plot_name, tolerance_dict, precursor_tolerances, frag_tolerances)
        all_matrices.append(z)
    base_matrix = np.zeroes(
        (
            len(precursor_tolerances),
            len(frag_tolerances),
        )
    )
    for m in all_matrices:
        pass

    # df = pd.DataFrame({'file': [os.path.basename(f) for f in files], 'offset': means})
    # df.to_csv('optimal_offsets.csv', index=False)


if __name__ == '__main__':
    main()
