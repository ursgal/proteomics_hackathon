#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import click
import pandas as pd
import ursgal
from scipy.stats import norm, normaltest
import sys


def draw_heatmap(tol_dict):
    pass


@click.command()
@click.argument('files', nargs=-1)
@click.argument('database')
@click.option('--force', '-f', is_flag=True, help="Print more output.")
def main(files, database, force):
    engine = 'msgfplus_v2019_04_18'
    engine = 'xtandem_alanine'
    engine = 'msfragger_20190222'
    params = {
        # 'use_pyqms_for_mz_calculation': True,
        'database': database,
        'scan_skip_modulo_step': 5,
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

    precursor_tolerances = [20, 15, 10, 5, 2, 1]
    frag_tolerances = [20, 15, 10, 5, 2, 1]
    for file in files:
        uc.params['machine_offset_in_ppm'] = 0
        uc.params['prefix'] = ''
        uc.params['precursor_mass_tolerance_unit'] = 'ppm'
        uc.params['precursor_mass_tolerance_plus'] = 300
        uc.params['precursor_mass_tolerance_minus'] = 300
        uc.params['frag_mass_tolerance'] = 700
        uc.params['frag_mass_tolerance_unit'] = 'ppm'
        # 'precursor_mass_tolerance_minus': 5,
        # 'precursor_mass_tolerance_plus': 5,

        # search = uc.search(file, engine=engine, force=True)

        raw_search_results = uc.execute_misc_engine(
            input_file = file,
            engine     = engine,
            force      = True,
        )
        csv_search_results = uc.convert(
            input_file = raw_search_results,
            engine = None,
            guess_engine = True,
            force      = force,
        )
        mapped_csv_search_results = uc.execute_misc_engine(
            input_file       = csv_search_results,
            engine           = uc.params['peptide_mapper_converter_version'],
            force            = force,
        )
        search = uc.execute_misc_engine(
            input_file       = mapped_csv_search_results,
            engine           = self.params['unify_csv_converter_version'],
            force            = force,
            merge_duplicates = True,
        )

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
                uc.params['prefix'] = f'corrected_mgf_{prec_tol}_{frag_tol}'
                uc.params['frag_mass_tolerance'] = frag_tol
                search = uc.search(corrected_mgf, engine=engine, force=force)
                validated = uc.validate(search, engine='percolator_3_4', force=force)
                # validated = uc.validate(search, engine='percolator_2_08', force=force)
                filtered = uc.execute_misc_engine(validated, engine='filter_csv_1_0_0', force=force)
                results = pd.read_csv(filtered)
                tolerance_dict[(prec_tol, frag_tol)] = len(results)
                from pprint import pprint
                pprint(tolerance_dict)

        draw_heatmap(dict)

    # df = pd.DataFrame({'file': [os.path.basename(f) for f in files], 'offset': means})
    # df.to_csv('optimal_offsets.csv', index=False)


if __name__ == '__main__':
    main()
