#!/usr/bin/env python3
import os
import ursgal
import sys
import glob
from collections import defaultdict as ddict
import csv


def main(database=None, folder=None):
    mzml_files = []
    for mzml in glob.glob(os.path.join(folder, "*.mzML")):
        # if 'TN_CSF_062617_03' not in mzml:
        #     continue
        mzml_files.append(mzml)

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
        "modifications": [],
        "-xmx": '32G',
        'psm_defining_colnames': [
            'Spectrum Title',
            'Sequence',
            'Modifications',
            'Mass Difference'
            'Charge',
            'Is decoy',
        ]
    }

    uc = ursgal.UController(params=params, profile="QExactive+", verbose=False)

    all_closed_search_files_nonfiltered = []
    all_closed_search_files_filtered = []
    for mzml in mzml_files:
        validated_results = os.path.join(
            folder,
            '{0}_pmap_unified_peptide_forest_1_0_0_validated.csv'.format(
                os.path.basename(mzml).replace('.mzML', '')
            )
        )
        uc.params['csv_filter_rules'] = [
            ['q-value_RF-reg','lte', 0.99],
        ]
        validated_results_sanitized = uc.execute_misc_engine(
            input_file = validated_results,
            engine='filter_csv',
        )
        all_closed_search_files_nonfiltered.append(validated_results_sanitized)

        uc.params['csv_filter_rules'] = [
            ['Is decoy', 'equals', 'false'],
            ['q-value_RF-reg','lte', 0.01],
            ['Conflicting uparam', 'contains_not', 'enzyme'],
        ]
        validated_results_filtered = uc.execute_misc_engine(
            input_file = validated_results_sanitized,
            engine='filter_csv',
        )
        all_closed_search_files_filtered.append(validated_results_filtered)

    merged_closed_search_nonfiltered = uc.execute_misc_engine(
        input_file = all_closed_search_files_nonfiltered,
        engine='merge_csvs',
        merge_duplicates=False,
        output_file_name='closed_search_all_results_sanitized.csv'
    )
    merged_closed_search_filtered = uc.execute_misc_engine(
        input_file = all_closed_search_files_filtered,
        engine='merge_csvs',
        merge_duplicates=False,
        output_file_name='closed_search_all_results_filtered.csv'
    )

    taggraph_nonfiltered = []
    taggraph_filtered = []
    for mzml in sorted(mzml_files)[:]:
        search_results = os.path.join(
            os.path.dirname(mzml),
            'tag_graph_1_8_0',
            'open_mod_{0}_tag_graph_1_8_0.csv'.format(
                os.path.basename(mzml).replace('.mzML', '')
            )
        )
        if os.path.exists(search_results) is False:
            search_results = os.path.join(
            os.path.dirname(mzml),
            'tag_graph_1_8_0',
            '{0}_tag_graph_1_8_0.csv'.format(
                os.path.basename(mzml).replace('.mzML', '')
            )
        )
        uc.params['prefix'] = 'open_mod'
        pmap_search_results = uc.execute_misc_engine(
            search_results,
            engine='upeptide_mapper_1_0_0',
        )
        unified_search_results = uc.execute_misc_engine(
            pmap_search_results,
            engine='unify_csv',
        )
        taggraph_nonfiltered.append(unified_search_results)

        uc.params['csv_filter_rules'] = [
            ['Is decoy', 'equals', 'false'],
            ['TagGraph:EM Probability','gte', 0.99],
            ['Conflicting uparam', 'contains_not', 'enzyme'],
        ]
        unified_search_results_filtered = uc.execute_misc_engine(
            input_file = unified_search_results,
            engine='filter_csv',
        )
        taggraph_filtered.append(unified_search_results_filtered)
    # exit()
    uc.params['prefix'] = ''
    merged_taggraph_nonfiltered = uc.execute_misc_engine(
        input_file = taggraph_nonfiltered,
        engine='merge_csvs',
        merge_duplicates=False,
        output_file_name='tag_graph_all_results.csv'
    )
    merged_taggraph_filtered = uc.execute_misc_engine(
        input_file = taggraph_filtered,
        engine='merge_csvs',
        merge_duplicates=False,
        output_file_name='tag_graph_all_results_accepted.csv'
    )

    open_search_nonfiltered = []
    open_search_filtered = []
    for mzml in mzml_files:
        ptminer_results = os.path.join(
            os.path.dirname(mzml),
            'ptminer_results',
            'open_mod_{0}___pmap_unified_peptide_forest_1_0_0_validated_accepted_ptminer_1_0.csv'.format(
                os.path.basename(mzml).replace('.mzML', '')
            )
        )
        open_search_nonfiltered.append(ptminer_results)

        uc.params['csv_filter_rules'] = [
            ['Is decoy', 'equals', 'false'],
            ['q-value_RF-reg','lte', 0.01],
            ['Conflicting uparam', 'contains_not', 'enzyme'],
        ]
        ptminer_results_filtered = uc.execute_misc_engine(
            input_file = ptminer_results,
            engine='filter_csv',
        )
        open_search_filtered.append(ptminer_results_filtered)

    merged_open_search_nonfiltered = uc.execute_misc_engine(
        input_file = open_search_nonfiltered,
        engine='merge_csvs',
        merge_duplicates=False,
        output_file_name='open_search_ptminer_results.csv'
    )
    merged_open_search_filtered = uc.execute_misc_engine(
        input_file = open_search_filtered,
        engine='merge_csvs',
        merge_duplicates=False,
        output_file_name='open_search_ptminer_results_accepted.csv'
    )

    pglyco_nonfiltered = []
    pglyco_filtered = []
    for mzml in mzml_files:
        unified_search_results = os.path.join(
            os.path.dirname(mzml),
            'pglyco_db_2_2_2',
            '{0}_pglyco_db_2_2_2_pglyco_fdr_2_2_2_pmap_unified.csv'.format(
                os.path.basename(mzml).replace('.mzML', '')
            )
        )
        pglyco_nonfiltered.append(unified_search_results)

        uc.params['csv_filter_rules'] = [
            ['Is decoy', 'equals', 'false'],
            ['q-value','lte', 0.01],
            ['Conflicting uparam', 'contains_not', 'enzyme'],
        ]
        unified_search_results_filtered = uc.execute_misc_engine(
            input_file = unified_search_results,
            engine='filter_csv',
        )
        pglyco_filtered.append(unified_search_results_filtered)

    merged_pglyco_nonfiltered = uc.execute_misc_engine(
        input_file = pglyco_nonfiltered,
        engine='merge_csvs',
        merge_duplicates=False,
        output_file_name='pglyco_all_results.csv'
    )

    merged_pglyco_filtered = uc.execute_misc_engine(
        input_file = pglyco_filtered,
        engine='merge_csvs',
        merge_duplicates=False,
        output_file_name='pglyco_all_results_accepted.csv'
    )






    merged_all_files_nonfiltered = uc.execute_misc_engine(
        input_file = [
            merged_closed_search_nonfiltered,
            merged_open_search_nonfiltered,
            merged_taggraph_nonfiltered,
            # merged_pglyco,
        ],
        engine='merge_csvs',
        merge_duplicates=False,
        output_file_name='closed_open_taggraph_merged_results.csv'
    )

    all_fieldnames = []
    spec_titles_with_ident = set()
    all_line_dicts = []
    with open(merged_all_files_nonfiltered , 'r') as file_object:
        csv_reader = csv.DictReader(file_object)
        fieldnames = csv_reader.fieldnames
        for fn in fieldnames:
            all_fieldnames.append(fn)
        for line_dict in csv_reader:
            spec_titles_with_ident.add(line_dict['Spectrum Title'])
            all_line_dicts.append(line_dict)

    with open(merged_pglyco_nonfiltered, 'r') as file_object:
        csv_reader = csv.DictReader(file_object)
        fieldnames = csv_reader.fieldnames
        for fn in fieldnames:
            if fn not in all_fieldnames:
                all_fieldnames.append(fn)
        for line_dict in csv_reader:
            if line_dict['Spectrum Title'] not in spec_titles_with_ident:
                spec_titles_with_ident.add(line_dict['Spectrum Title'])
                all_line_dicts.append(line_dict)

    csv_kwargs = {}
    if sys.platform == 'win32':
        csv_kwargs['lineterminator'] = '\n'
    else:
        csv_kwargs['lineterminator'] = '\r\n'

    out_file_path = os.path.join(
        os.path.dirname(merged_all_files_nonfiltered),
        'All_merged_results_unfiltered.csv'
    )
    with open(out_file_path, 'w') as out_csv:
        csv_writer = csv.DictWriter(
            out_csv,
            all_fieldnames,
            **csv_kwargs
        )
        csv_writer.writeheader()
        for out_dict in all_line_dicts:
            csv_writer.writerow(out_dict)



    merged_all_files_filtered = uc.execute_misc_engine(
        input_file = [
            merged_closed_search_filtered,
            merged_open_search_filtered,
            merged_taggraph_filtered,
            # merged_pglyco,
        ],
        engine='merge_csvs',
        merge_duplicates=False,
        output_file_name='closed_open_taggraph_merged_results_accepted.csv'
    )

    all_fieldnames = []
    spec_titles_with_ident = set()
    all_line_dicts = []
    with open(merged_all_files_filtered , 'r') as file_object:
        csv_reader = csv.DictReader(file_object)
        fieldnames = csv_reader.fieldnames
        for fn in fieldnames:
            all_fieldnames.append(fn)
        for line_dict in csv_reader:
            spec_titles_with_ident.add(line_dict['Spectrum Title'])
            all_line_dicts.append(line_dict)

    with open(merged_pglyco_filtered, 'r') as file_object:
        csv_reader = csv.DictReader(file_object)
        fieldnames = csv_reader.fieldnames
        for fn in fieldnames:
            if fn not in all_fieldnames:
                all_fieldnames.append(fn)
        for line_dict in csv_reader:
            if line_dict['Spectrum Title'] not in spec_titles_with_ident:
                spec_titles_with_ident.add(line_dict['Spectrum Title'])
                all_line_dicts.append(line_dict)

    csv_kwargs = {}
    if sys.platform == 'win32':
        csv_kwargs['lineterminator'] = '\n'
    else:
        csv_kwargs['lineterminator'] = '\r\n'

    out_file_path = os.path.join(
        os.path.dirname(merged_all_files_nonfiltered),
        'All_merged_results_accepted.csv'
    )
    with open(out_file_path, 'w') as out_csv:
        csv_writer = csv.DictWriter(
            out_csv,
            all_fieldnames,
            **csv_kwargs
        )
        csv_writer.writeheader()
        for out_dict in all_line_dicts:
            csv_writer.writerow(out_dict)

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
    )