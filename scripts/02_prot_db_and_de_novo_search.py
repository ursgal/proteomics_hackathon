#!/usr/bin/env python3
import os
import pandas as pd
import ursgal
import sys
import glob
from collections import defaultdict as ddict
import csv


def main(folder=None, database=None, offset_file=None, approach="prot_db"):
    mzML_files = []
    for mzml in glob.glob(os.path.join(folder, "*.mzML")):
        mzML_files.append(mzml)

    params = {
        "use_pyqms_for_mz_calculation": True,
        "cpus": 8,
        "database": database,  # TREMBL + SwissProt
        "precursor_mass_tolerance_unit": "ppm",
        "precursor_mass_tolerance_plus": 20,
        "precursor_mass_tolerance_minus": 20,
        "frag_mass_tolerance_unit": "ppm",
        "frag_mass_tolerance": 0.005,
        "frag_mass_tolerance_unit": "da",
        "csv_filter_rules": [
            ["Is decoy", "equals", "false"],
            ["PEP", "lte", 0.01],
        ],
        "modifications": [
            "M,opt,any,Oxidation",
            "*,opt,Prot-N-term,Acetyl",
            'N,opt,any,Deamidated',
            'Q,opt,any,Deamidated',
            "C,fix,any,Carbamidomethyl",
        ],
        "peptide_mapper_class_version": "UPeptideMapper_v4",
        "-xmx": '32G'
    }

    uc = ursgal.UController(params=params, profile="QExactive+", verbose=False)

    offset_dict = {}
    with open(offset_file, "r") as offset_input:
        csv_reader = csv.DictReader(offset_input)
        for line_dict in csv_reader:
            offset_dict[line_dict["File name"]] = float(line_dict["Optimum Precursor offset"])

    for mzml in sorted(mzML_files)[:]:
        if 'TN_CSF_062617_01.mzML' in mzml:
            continue
        percolator_results = []
        unvalidated_results = []
        mzml_basename = os.path.basename(mzml)
        uc.params['machine_offset_in_ppm'] = offset_dict[mzml_basename]
        if approach == "prot_db":
            for prot_db_engine in [
                # 'msfragger_2_3',
                # 'msgfplus_v2019_07_03',
                # 'mascot_2_6_0',
                # "xtandem_alanine",
                # "omssa_2_1_9",
                "msamanda_2_0_0_14665",
            ]:
                unified_search_results = uc.search(mzml, engine=prot_db_engine,)
                unvalidated_results.append(unified_search_results)

                vd_engine = "percolator_3_4_0"
                validated_results = uc.validate(
                    unified_search_results, engine=vd_engine,
                )
                percolator_results.append(validated_results)

        if approach == "de_novo":
            de_novo_results = []
            for de_novo_engine in [
                # 'pnovo_3_1_3',
                'novor_1_05',
                'deepnovo_v2',
            ]:
                # raw_file = mzml.replace('.mzML', '.raw')
                # xtracted_file = uc.convert(
                #     input_file = raw_file,
                #     engine = 'pparse_2_2_1',
                #     # force = True,
                # )
                # search_result = uc.search_mgf(
                #     input_file = xtracted_file,
                #     engine = de_novo_engine,
                # )
                mgf_file = uc.convert(
                    input_file = mzml,
                    engine = 'mzml2mgf_2_0_0',
                    # force = True,
                )
                # unified_search_results = uc.execute_misc_engine(
                #     input_file = search_result,
                #     engine='unify_csv'
                # )
                # uc.params['precursor_max_mass'] = 1000
                if de_novo_engine == 'novor_1_05':
                    uc.params["modifications"] = [
                        "M,opt,any,Oxidation",
                        "*,opt,Prot-N-term,Acetyl",
                        "C,fix,any,Carbamidomethyl",
                    ]
                else:
                    uc.params["modifications"] = [
                        "M,opt,any,Oxidation",
                        "*,opt,Prot-N-term,Acetyl",
                        'N,opt,any,Deamidated',
                        'Q,opt,any,Deamidated',
                        "C,fix,any,Carbamidomethyl",
                    ]
                uc.params['deepnovo_use_lstm'] = False
                search_result = uc.search_mgf(
                    input_file = mgf_file,
                    engine = de_novo_engine
                )
                unified_search_results = uc.execute_misc_engine(
                    input_file = search_result,
                    engine='unify_csv'
                )
                de_novo_results.append(unified_search_results)

            uc.params['prefix'] = 'de_novo'
            merged_all_perc = uc.execute_misc_engine(
                input_file = de_novo_results,
                engine='merge_csvs',
                merge_duplicates=False,
            )
            uc.params['prefix'] = ''


    #     vd_engine = 'peptide_forest'
    #     peptide_forest_alone = uc.valdiate(vd_engine, unvalidated_results)
    #     peptide_forest_alone = uc.valdiate(vd_engine, de_novo_results)
    #     peptide_forest_percolator = uc.valdiate(vd_engine, percolator_results)

    #     venn_diagram = uc.visualize([percolator_results, peptide_forest_alone, peptide_forest_percolator])

    # params = {optimized_from_sweep, wide_search_arams}
    # for file in files:
    #     percolator_results = []
    #     unvalidated_results = []
    #     for open_mod_engine in [
    #         'msfragger',
    #         'moda',
    #         'pipi',
    #         'taggraph',
    #     ]:
    #         unified_search_results = uc.search(file, prot_db_engine)
    #         unvalidated_results.append(unified_search_results)

    #         vd_engine = 'percolator_3_4'
    #         validated_results = uc.valdiate(vd_engine, unified_search_results)
    #         percolator_results.append(validated_results)

    #     vd_engine = 'peptide_forest'
    #     peptide_forest_alone = uc.valdiate(vd_engine, unvalidated_results)
    #     peptide_forest_percolator = uc.valdiate(vd_engine, percolator_results)

    #     venn_diagram = uc.visualize([percolator_results, peptide_forest_alone, peptide_forest_percolator])


if __name__ == '__main__':
    main(
        folder=sys.argv[1],
        database=sys.argv[2],
        offset_file=sys.argv[3],
        approach=sys.argv[4],
    )


