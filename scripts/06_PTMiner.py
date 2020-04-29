#!/usr/bin/env python3
import os
import ursgal
import sys
import glob
from collections import defaultdict as ddict
import csv


def main(folder=None, database=None):
    if folder.endswith('.mgf'):
        mgf_files = [folder]
    else:    
        mgf_files = []
        for mzml in glob.glob(os.path.join(folder, "*.mgf")):
            mgf_files.append(mzml)

    params = {
        "use_pyqms_for_mz_calculation": True,
        'database': database,
        "cpus": 8,
        "precursor_mass_tolerance_unit": "ppm",
        "precursor_mass_tolerance_plus": 20,
        "precursor_mass_tolerance_minus": 20,
        "frag_mass_tolerance_unit": "ppm",
        "frag_mass_tolerance": 0.005,
        "frag_mass_tolerance_unit": "da",
        "modifications": [
            "C,fix,any,Carbamidomethyl",
        ],
        "-xmx": '32G',
        'psm_defining_colnames': [
            'Spectrum Title',
            'Sequence',
            'Modifications',
            'Charge',
            'Is decoy',
        ]
    }

    uc = ursgal.UController(params=params, profile="QExactive+", verbose=False)

    # all_mgf_files = []
    # for mzml in sorted(mzML_files)[:]:
    #     all_mgf_files.append(mzml.replace('.mzML', '.mgf'))
    #     mzml_basename = os.path.basename(mzml)

    all_ptminer_unfiltered = []
    all_ptminer_filtered = []
    all_pf_filtered = []
    for mgf_file in sorted(mgf_files)[:]:
        peptife_forest_results = os.path.join(
            os.path.dirname(mgf_file),
            os.path.basename(mgf_file).replace(
                '.mgf', '___pmap_unified_peptide_forest_1_0_0_validated.csv'
            ),
        )

        uc.params["csv_filter_rules"] = [
            ["q-value_RF-reg", "lte", 0.99],
            ["Conflicting uparam", "contains_not", "enzyme"],
        ]
        filtered_peptide_forest = uc.execute_misc_engine(
            input_file=peptife_forest_results,
            engine='filter_csv_1_0_0',
        )
        # all_pf_filtered.append(filtered_peptide_forest)

        uc.params.update({
            'mgf_input_files_list' : [mgf_file],
            'validation_score_field': 'q-value_RF-reg',
            'bigger_scores_better': False,
        })
        ptminer_peptide_forest = uc.validate(
            input_file=filtered_peptide_forest,
            engine='ptminer',
            # force=True,
        )
        all_ptminer_unfiltered.append(ptminer_peptide_forest)

        uc.params["csv_filter_rules"] = [
            ["Is decoy", "equals", "false"],
            ["q-value_RF-reg", "lte", 0.01],
            ["Conflicting uparam", "contains_not", "enzyme"],
        ]

        filtered_ptminer_peptide_forest = uc.execute_misc_engine(
            input_file=ptminer_peptide_forest,
            engine='filter_csv_1_0_0',
        )
        all_ptminer_filtered.append(filtered_ptminer_peptide_forest)

    merged_all_pf_filtered = uc.execute_misc_engine(
        input_file = all_ptminer_filtered,
        engine='merge_csvs',
        merge_duplicates=False,
    )
    merged_all_pf_unfiltered = uc.execute_misc_engine(
        input_file = all_ptminer_unfiltered,
        engine='merge_csvs',
        merge_duplicates=False,
    )

    # uc.params.update({
    #     'mgf_input_files_list' : all_mgf_files,
    #     'validation_score_field': 'combined PEP',
    #     'bigger_scores_better': False,
    # })

    # ptminer_all_combined_pep = uc.validate(
    #     input_file=merged_csv_combined_pep,
    #     engine='ptminer',
    #     # force=True,
    # )

    # uc.params['csv_filter_rules'] = [
    #     ['Is decoy', 'equals', 'false'],
    #     ['combined PEP','lte', 0.01],
    #     ['Conflicting uparam', 'contains_not', 'enzyme'],
    # ]

    # filtered_ptminer_all_combined_pep = uc.execute_misc_engine(
    #     input_file=ptminer_all_combined_pep,
    #     engine='filter_csv_1_0_0',
    # )

    # uc.params.update({
    #     'visualization_column_names': ['Modifications', 'Sequence'],
    #     'visualization_label_positions':{
    #         '0': 'Percolator',
    #         '1': 'Combined PEP',
    #     }
    # })
    # uc.visualize(
    #     input_files=[
    #         filtered_ptminer_all_perc,
    #         filtered_ptminer_all_combined_pep,
    #     ],
    #     engine='venndiagram_1_1_0',
    #     output_file_name='open_search_ptminer_venndiagram_peptides.svg'
    # )

    # uc.params['visualization_column_names'] = ['Spectrum Title', 'Sequence', 'Modifications']
    # uc.visualize(
    #     input_files=[
    #         filtered_ptminer_all_perc,
    #         filtered_ptminer_all_combined_pep,
    #     ],
    #     engine='venndiagram_1_1_0',
    #     output_file_name='open_search_ptminer_venndiagram_psms.svg'
    # )

if __name__ == '__main__':
    main(
        folder=sys.argv[1],
        database=sys.argv[2],
        # merged_csv_perc=sys.argv[3],
        # merged_csv_combined_pep=sys.argv[4]
    )