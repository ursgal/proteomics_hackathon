#!/usr/bin/env python3
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import click
import pandas as pd
import ursgal
from scipy.stats import norm, normaltest
import sys


def draw_heatmap(plot_name, prec_tolerances, frag_tolerances, tol_dict=None, z=None,):
    # breakpoint()
    if tol_dict is not None and z is None:
        z = np.zeros(
            (len(prec_tolerances), len(frag_tolerances))
        )
        for i, prec_tol in enumerate(prec_tolerances):
            for j, frag_tol in enumerate(frag_tolerances):
                z[i, j] = tol_dict[(prec_tol, frag_tol)]
    # z /= z.max()
    if z is not None:
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
    engine = 'xtandem_alanine'
    engine = 'msfragger_2_3'
    params = {
        # 'use_pyqms_for_mz_calculation': True,
        'cpus': 6,
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
    frag_tolerances      = [20  , 17.5, 15  , 12.5, 10, 5]
    all_matrices = []
    csv_lines= []
    for file in files:
        # if ('_01.mzML' in file) or ('_47.mzML' in file):# or ('_11.mzML' in file) or ('_19.mzML' in file):
        #     continue
        uc.params['machine_offset_in_ppm'] = 0
        uc.params['prefix'] = ''
        uc.params['precursor_mass_tolerance_unit'] = 'ppm'
        uc.params['precursor_mass_tolerance_plus'] = 300
        uc.params['precursor_mass_tolerance_minus'] = 300
        uc.params['frag_mass_tolerance'] = 700
        uc.params['frag_mass_tolerance_unit'] = 'ppm'
        # 'precursor_mass_tolerance_minus': 5,
        # 'precursor_mass_tolerance_plus': 5,
        try:
            search = uc.search(file, engine=engine, force=force)
        except:
            continue
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

        max_res = 0
        max_param = None
        for prec_tol in precursor_tolerances:
            uc.params['precursor_mass_tolerance_plus'] = prec_tol
            uc.params['precursor_mass_tolerance_minus'] = prec_tol
            for frag_tol in frag_tolerances:
                # if prec tol > frag_tol:
                #     continue
                uc.params['prefix'] = f'corrected_mgf_{prec_tol}_{frag_tol}'
                uc.params['frag_mass_tolerance'] = frag_tol
                # try:
                search = uc.search(corrected_mgf, engine=engine, force=force)
                validated = uc.validate(search, engine='percolator_3_4', force=force)
                # validated = uc.validate(search, engine='percolator_2_08', force=force)
                filtered = uc.execute_misc_engine(validated, engine='filter_csv_1_0_0', force=force)
                results = pd.read_csv(filtered)
                # except:
                    # results = []
                tolerance_dict[(prec_tol, frag_tol)] = len(results)
                if len(results) > max_res:
                    max_res = len(results)
                    max_param = (prec_tol, frag_tol)
                file_basename = os.path.splitext(os.path.basename(file))[0]
                plot_name = f'{file_basename}_{engine}_heatmap.png'
        csv_lines.append(
            {
                'File name': os.path.basename(file),
                'Optimum Precursor offset': max_param[0],
                'Optimum Precursor unit': uc.params['precursor_mass_tolerance_unit'],
                'Optimal Fragment offset': max_param[1],
                'Optimal Fragment unit': uc.params['frag_mass_tolerance_unit'],
            }
        )
        print(tolerance_dict)
        z = draw_heatmap(plot_name, precursor_tolerances, frag_tolerances, tol_dict=tolerance_dict)
        all_matrices.append(z)
    field_names = [
        'File name',
        'Optimum Precursor offset',
        'Optimum Precursor unit',
        'Optimal Fragment offset',
        'Optimal Fragment unit',
    ]
    with open('tolerances.csv', 'wt') as fout:
        writer = csv.DictWriter(fout, fieldnames=field_names)
        writer.writeheader()
        for line in csv_lines:
            writer.writerow(line)

    combined_matrix = sum(all_matrices)
    draw_heatmap(f'{engine}_heatmap_all_files.png', precursor_tolerances, frag_tolerances, z=combined_matrix)

    # df = pd.DataFrame({'file': [os.path.basename(f) for f in files], 'offset': means})
    # df.to_csv('optimal_offsets.csv', index=False)


if __name__ == '__main__':
    main()
