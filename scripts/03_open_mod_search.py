#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os
import csv


def main(folder=None, target_decoy_database=None, offset_file=None):
    mzML_files = []
    for mzml in glob.glob(os.path.join('{0}'.format(folder), '*.mzML')):
        mzML_files.append(mzml)

    mass_spectrometer = 'QExactive+'

    # We specify all search engines and validation engines that we want to use in a list
    # (version numbers might differ on windows or mac):
    search_engines = [
        'msfragger_2_3',
        'moda_v1_61',
        'pipi_1_4_6',
    ]

    validation_engines = [
        'percolator_3_4_0',
    ]

    # Initializing the Ursgal UController class with
    # our specified mass spectrometer

    params = {
        "use_pyqms_for_mz_calculation": True,
        'database': target_decoy_database,
        'modifications' : ['C,fix,any,Carbamidomethyl'],
        'csv_filter_rules': [
            ['Is decoy', 'equals', 'false'],
            ['PEP', 'lte', 0.01],
        ],
        'enzyme' : 'trypsin',
        'frag_mass_tolerance_unit' : 'da',
        'frag_mass_tolerance' : 0.005,
        'precursor_mass_tolerance_unit' : 'ppm',
        'precursor_mass_tolerance_plus' : 5,
        'precursor_mass_tolerance_minus' : 5,
        'moda_high_res' : False,
        'max_mod_size' : 2000,
        'min_mod_size' : -230,
        'precursor_true_units' : 'ppm',
        'precursor_true_tolerance' : 5,
        'remove_temporary_files' : False,
        "prefix": "open_mod",
    }

    uc = ursgal.UController(
        profile="QExactive+",
        params=params,
    )

    offset_dict = {}
    with open(offset_file, "r") as offset_input:
        csv_reader = csv.DictReader(offset_input)
        for line_dict in csv_reader:
            offset_dict[line_dict["File name"]] = float(line_dict["Optimum Precursor offset"])

    # complete workflow:
    # every spectrum file is searched with every search engine,
    # results are validated and filtrated for targets and PEP <= 0.01 (for each engine seperately).
    # In the end, all filtered results from all spectrum files are merged
    for validation_engine in validation_engines:
        all_engines_results = []
        for search_engine in search_engines :
            result_all_files_each_engine = []
            for spec_file in sorted(mzML_files)[:]:
                if 'TN_CSF_062617_01.mzML' in spec_file:
                    continue
                mzml_basename = os.path.basename(spec_file)
                uc.params['machine_offset_in_ppm'] = offset_dict[mzml_basename]
                #1. convert to MGF
                mgf_file = uc.convert(
                    input_file=spec_file,
                    engine = 'mzml2mgf_2_0_0',
                    # force =True,
                    )
                # continue
                if search_engine == 'msfragger_2_3':
                    uc.params.update(
                        {'precursor_mass_tolerance_unit' : 'da',
                        'precursor_mass_tolerance_plus' : 2000,
                        'precursor_mass_tolerance_minus' : 230}
                    )
                # if search_engine == 'pipi_1_4_5':
                #     uc.params.update({'database':target_decoy_database,})

                #2. do the actual search
                raw_search_results=uc.search_mgf(
                    input_file = mgf_file,
                    engine = search_engine,
                    )

                uc.params.update(
                    {
                        'database': target_decoy_database,
                        'precursor_mass_tolerance_unit': 'ppm',
                        'precursor_mass_tolerance_plus': 5,
                        'precursor_mass_tolerance_minus': 5
                    }
                )

                #3. convert files to csv
                csv_search_results= uc.convert(
                    input_file=raw_search_results,
                    engine = None,
                    guess_engine = True,
                    )

                #4. protein mapping.
                mapped_csv_search_results = uc.execute_misc_engine(
                    input_file       = csv_search_results,
                    #output_file_name = output_file_name,
                    engine           = 'upeptide_mapper_1_0_0',
                )

                # 5. Convert csv to unified ursgal csv format:
                unified_search_results = uc.execute_misc_engine(
                    input_file       = mapped_csv_search_results,
                    #output_file_name = output_file_name,
                    engine           = 'unify_csv_1_0_0',
                    merge_duplicates = False,
                )

                validated_csv = uc.validate(
                    input_file=unified_search_results,
                    engine=validation_engine,
                    )

                #validated_results.append(validated_csv)

                filtered_validated_results = uc.execute_misc_engine(
                    input_file=validated_csv,
                    engine='filter_csv_1_0_0',
                )
                result_all_files_each_engine.append(filtered_validated_results)

            merged_files_from_each_engine = uc.execute_misc_engine(
                input_file=result_all_files_each_engine,
                engine='merge_csvs_1_0_0',
                merge_duplicates=True,
                )

            all_engines_results.append(merged_files_from_each_engine)

        merged_files_from_all_engine = uc.execute_misc_engine(
            input_file=all_engines_results,
            engine='merge_csvs_1_0_0',
            merge_duplicates=True,
            )


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(main.__doc__)
        sys.exit(1)
    main(
        folder=sys.argv[1],
        target_decoy_database=sys.argv[2],
        offset_file=sys.argv[3],
    )
