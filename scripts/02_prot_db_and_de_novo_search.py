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
        # 'cpus': 8,
        "database": database,  # TREMBL + SwissProt
        "precursor_mass_tolerance_unit": "ppm",
        "precursor_mass_tolerance_plus": 20,
        "precursor_mass_tolerance_minus": 20,
        "frag_mass_tolerance_unit": "ppm",
        "frag_mass_tolerance": 0.005,
        "frag_mass_tolerance_unit": "da",
        "csv_filter_rules": [
            ["Is decoy", "equals", "false"],
            ["q-value", "lte", 0.01],
        ],
        "modifications": [
            "M,opt,any,Oxidation",
            "*,opt,Prot-N-term,Acetyl",
            "C,fix,any,Carbamidomethyl",
        ],
        "peptide_mapper_class_version": "UPeptideMapper_v4",
    }

    uc = ursgal.UController(params=params, profile="QExactive+", verbose=False)

    offset_dict = {}
    with open(offset_file, "r") as offset_input:
        csv_reader = csv.DictReader(offset_input)
        for line_dict in csv_reader:
            offset_dict[line_dict["File name"]] = line_dict["Optimum Precursor offset"]

    for mzml in mzML_files:
        percolator_results = []
        unvalidated_results = []
        uc.params["machine_offset_in_ppm"] = offset_dict[mzml]
        if approach == "prot_db":
            for prot_db_engine in [
                'msgfplus_v2019_07_03',
                'msfragger_2_3',
                'mascot_2_6_0',
                "xtandem_alanine",
                "omssa_2_1_9",
                "msamanda_2_0_0_13723",
            ]:
                unified_search_results = uc.search(mzml, engine=prot_db_engine,)
                unvalidated_results.append(unified_search_results)

                vd_engine = "percolator_3_4_0"
                validated_results = uc.valdiate(
                    unified_search_results, engine=vd_engine,
                )
                percolator_results.append(validated_results)

        if approach == "de_novo":
            de_novo_results = []
            for de_novo_engine in [
                'pnovo_3_1_3',
                'novor_1_05',
                # Do deepnovo
                'deepnovo_0_0_1'
            ]:
                unified_search_results = uc.search(file, engine=de_novo_engine)
                de_novo_results.append(unified_search_results)

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


if __name__ == "__main__":
    main(
        folder=sys.argv[1],
        database=sys.argv[2],
        offset_file=sys.argv[3],
        approach=sys.argv[4],
    )
