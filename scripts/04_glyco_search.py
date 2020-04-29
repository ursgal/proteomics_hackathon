#!/usr/bin/env python3.4
# encoding: utf-8

import ursgal
import os
import sys
import shutil
import csv
import glob


def main(folder=None, database=None, offset_file=None):
    # mzML_files = []
    # for mzml in glob.glob(os.path.join(folder, "*.mzML")):
    #     mzML_files.append(mzml)

    uc = ursgal.UController(
        profile='QExactive+',
        params={
            'database': database,
            'modifications': [
                'M,opt,any,Oxidation',        # Met oxidation
                'C,fix,any,Carbamidomethyl',  # Carbamidomethylation
                '*,opt,Prot-N-term,Acetyl'    # N-Acteylation
            ],
            'peptide_mapper_class_version': 'UPeptideMapper_v4',
            'enzyme': 'trypsin',
            'frag_mass_tolerance'       : 20,
            'frag_mass_tolerance_unit'  : 'ppm',
            'precursor_mass_tolerance_plus' : 5,
            'precursor_mass_tolerance_minus' : 5,
            'aa_exception_dict' : {},
            "use_pyqms_for_mz_calculation": True,
            "cpus": 8,
            "precursor_mass_tolerance_unit": "ppm",
            "precursor_mass_tolerance_plus": 20,
            "precursor_mass_tolerance_minus": 20,
            "frag_mass_tolerance_unit": "ppm",
            "frag_mass_tolerance": 0.005,
            "frag_mass_tolerance_unit": "da",
            "csv_filter_rules": [
                ["Is decoy", "equals", "false"],
                ["q-value", "lte", 0.01],
                ["Conflicting uparam", "contains_not", "enzyme"]
            ],
            "-xmx": '32G'
        }
    )

    offset_dict = {}
    with open(offset_file, "r") as offset_input:
        csv_reader = csv.DictReader(offset_input)
        for line_dict in csv_reader:
            offset_dict[line_dict["File name"]] = float(line_dict["Optimum Precursor offset"])

    raw_files = []
    for raw in glob.glob(os.path.join('{0}'.format(folder), '*.raw')):
        raw_files.append(raw)

    filtered_results_list = []
    for raw_file in raw_files:
        # xtracted_file = uc.convert(
        #     input_file = raw_file,
        #     engine = 'pparse_2_2_1',
        #     # force = True,
        # )

        # mzml_basename = os.path.basename(raw_file.replace('.raw', '.mzML'))
        # uc.params['machine_offset_in_ppm'] = offset_dict[mzml_basename]
        # mgf_file = uc.convert(
        #     input_file = raw_file.replace('.raw', '.mzML'),
        #     engine = 'mzml2mgf_2_0_0',
        #     # force = True,
        # )

        # search_result = uc.search_mgf(
        #     input_file = xtracted_file,
        #     engine = 'pglyco_db_2_2_2',
        # )
            
        # mapped_results = uc.execute_misc_engine(
        #     input_file=search_result,
        #     engine='upeptide_mapper',
        # )

        # unified_search_results = uc.execute_misc_engine(
        #     input_file = mapped_results,
        #     engine='unify_csv'
        # )

        # validated_file = uc.validate(
        #     input_file=search_result,
        #     engine='pglyco_fdr_2_2_2',
        # )

        # mapped_validated_results = uc.execute_misc_engine(
        #     input_file=validated_file,
        #     engine='upeptide_mapper',
        # )

        # unified_validated_results = uc.execute_misc_engine(
        #     input_file = mapped_validated_results,
        #     engine='unify_csv',
        # )

        mzml_basename = os.path.basename(raw_file.replace('.raw', '.mzML'))
        unified_validated_results = os.path.join(
            os.path.dirname(raw_file),
            'pglyco_db_2_2_2',
            mzml_basename.replace('.mzML', '_pglyco_db_2_2_2_pglyco_fdr_2_2_2_pmap_unified.csv')
        )

        filtered_validated_results = uc.execute_misc_engine(
            input_file = unified_validated_results,
            engine='filter_csv',
        )

        filtered_results_list.append(filtered_validated_results)

    results_all_files = uc.execute_misc_engine(
        input_file = filtered_results_list,
        engine='merge_csv',
    )

    return


if __name__ == '__main__':
    main(
        folder=sys.argv[1],
        database=sys.argv[2],
        offset_file=sys.argv[3],
    )
