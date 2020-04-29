#!/usr/bin/env python3
import os
import ursgal
import sys
import glob
from collections import defaultdict as ddict
import csv


def main(folder=None, database=None):
    mzML_files = []
    for mzml in glob.glob(os.path.join(folder, "*.mzML")):
        mzML_files.append(mzml)

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
            "M,opt,any,Oxidation",
            "*,opt,Prot-N-term,Acetyl",
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

    all_percolator_results = []
    all_combined_pep_results = []
    all_peptide_forest_alone = []
    all_peptide_forest_percolator = []
    all_peptide_forest_de_novo = []
    all_open_validated = []
    all_open_combined_pep = []
    all_open_peptide_forest_alone = []
    all_open_peptide_forest_percolator = []
    all_mgf_files = []
    for mzml in sorted(mzML_files)[:]:
        all_mgf_files.append(mzml.replace('.mzML', '.mgf'))
        mzml_basename = os.path.basename(mzml)
        ##################################
        # Prot DB section
        ##################################

        percolator_results = []
        unvalidated_results = []
        for prot_db_engine in [
            'msfragger_2_3',
            'msgfplus_v2019_07_03',
            # 'mascot_2_6_0',
            "xtandem_alanine",
            "omssa_2_1_9",
            "msamanda_2_0_0_14665",
        ]:
            if prot_db_engine in [
                'omssa_2_1_9',
                "xtandem_alanine",
            ]:
                prefix = 'main_search_'
                uc.params['prefix'] = prefix
            else:
                prefix = ''
                uc.params['prefix'] = prefix
            pmap_search_results = os.path.join(
                os.path.dirname(mzml),
                prot_db_engine,
                prefix + mzml_basename.replace(
                    '.mzML', '_{0}_pmap.csv'.format(prot_db_engine)
                )
            )
            if os.path.exists(pmap_search_results) is False:
                search_results = os.path.join(
                    os.path.dirname(mzml),
                    prot_db_engine,
                    prefix + mzml_basename.replace(
                        '.mzML', '_{0}.csv'.format(prot_db_engine)
                    )
                )
                pmap_search_results = uc.execute_misc_engine(
                    search_results,
                    engine='upeptide_mapper_1_0_0',
                )    
            unified_search_results = uc.execute_misc_engine(
                pmap_search_results,
                engine='unify_csv',
            )

            unvalidated_results.append(unified_search_results)
            
            validated_results = uc.validate(
                unified_search_results,
                engine='percolator_3_4',
            )
            percolator_results.append(validated_results)

            # uc.params["csv_filter_rules"] = [
            #     ["Is decoy", "equals", "false"],
            #     ["PEP", "lte", 0.01],
            #     ["Conflicting uparam", "contains_not", "enzyme"],
            # ]
            # filtered_validated_results = uc.execute_misc_engine(
            #     input_file = validated_results,
            #     engine='filter_csv',
            # )
            # all_percolator_results.append(filtered_validated_results)

        # uc.params['prefix'] = mzml_basename.replace('.mzML', '')
        # combined_pep = uc.combine_search_results(
        #     input_files=percolator_results,
        #     engine='combine_pep_1_0_0',
        # )
        # uc.params['csv_filter_rules'] = [
        #     ['Is decoy', 'equals', 'false'],
        #     ['combined PEP','lte', 0.01],
        #     ['Conflicting uparam', 'contains_not', 'enzyme'],
        # ]
        # filtered_combined_results = uc.execute_misc_engine(
        #     input_file = combined_pep,
        #     engine='filter_csv',
        # )
        # all_combined_pep_results.append(filtered_combined_results)
        # uc.params['prefix'] = ''

        # ##################################
        # # De Novo section
        # ##################################

        # de_novo_results = []
        # for de_novo_engine in [
        #     # 'pnovo_3_1_3',
        #     'novor_1_05',
        #     'deepnovo_v2'
        # ]:
        #     unified_search_results = os.path.join(
        #         os.path.dirname(mzml),
        #         de_novo_engine,
        #         mzml_basename.replace('.mzML', '_{0}_pmap_unified.csv'.format(de_novo_engine))
        #     )
        #     de_novo_results.append(unified_search_results)

        # ##################################
        # # Open Mod section
        # ##################################

        # open_validated = []
        # open_unvalidated = []
        # for open_mod_engine in [
        #     'msfragger_2_3',
        #     'moda_v1_61',
        #     'pipi_1_4_6',
        #     # 'taggraph',
        # ]:

        #     unified_search_results = os.path.join(
        #         os.path.dirname(mzml),
        #         open_mod_engine,
        #         'open_mod_' + mzml_basename.replace(
        #             '.mzML', '_{0}_pmap_unified.csv'.format(open_mod_engine)
        #         )
        #     )
        #     open_unvalidated.append(unified_search_results)

        #     validated_results = uc.validate(
        #         unified_search_results,
        #         engine='percolator_3_4',
        #     )
        #     open_validated.append(validated_results)

        #     uc.params["csv_filter_rules"] = [
        #         ["Is decoy", "equals", "false"],
        #         ["PEP", "lte", 0.01],
        #         ["Conflicting uparam", "contains_not", "enzyme"],
        #     ]
        #     filtered_validated_results = uc.execute_misc_engine(
        #         input_file = validated_results,
        #         engine='filter_csv',
        #     )
        #     all_open_validated.append(filtered_validated_results)

        # uc.params['prefix'] = 'open_mod_' + mzml_basename.replace('.mzML', '')
        # open_combined_pep = uc.combine_search_results(
        #     input_files=open_validated,
        #     engine='combine_pep_1_0_0',
        # )
        # uc.params['csv_filter_rules'] = [
        #     ['Is decoy', 'equals', 'false'],
        #     ['combined PEP','lte', 0.01],
        #     ['Conflicting uparam', 'contains_not', 'enzyme'],
        # ]
        # open_filtered_combined_pep = uc.execute_misc_engine(
        #     input_file = open_combined_pep,
        #     engine='filter_csv',
        # )
        # all_open_combined_pep.append(open_filtered_combined_pep)
        # uc.params['prefix'] = ''


        ##################################
        # Peptide Forest section
        ##################################

        uc.params['peptide_forest_initial_engine'] = 'msgfplus_v2019_07_03'
        uc.params['peptide_forest_file_params'] = {}
        uc.params['prefix'] = mzml_basename.replace('.mzML', '')

        validated_peptide_forest_alone = uc.validate(
            input_file=unvalidated_results,
            engine='peptide_forest',
        )
        uc.params['csv_filter_rules'] = [
            ['Is decoy', 'equals', 'false'],
            ['q-value_RF-reg','lte', 0.01],
            ['Conflicting uparam', 'contains_not', 'enzyme'],
        ]
        filtered_peptide_forest_alone = uc.execute_misc_engine(
            input_file = validated_peptide_forest_alone,
            engine='filter_csv',
        )
        all_peptide_forest_alone.append(filtered_peptide_forest_alone)
        uc.params['prefix'] = ''
        # uc.params['prefix'] = 'pf_perc_' + mzml_basename.replace('.mzML', '')
        # validated_peptide_forest_percolator = uc.validate(
        #     input_file=percolator_results,
        #     engine='peptide_forest',
        # )
        # uc.params['csv_filter_rules'] = [
        #     ['Is decoy', 'equals', 'false'],
        #     ['q-value_RF-reg','lte', 0.01],
        #     ['Conflicting uparam', 'contains_not', 'enzyme'],
        # ]
        # filtered_peptide_forest_perc = uc.execute_misc_engine(
        #     input_file = validated_peptide_forest_percolator,
        #     engine='filter_csv',
        # )
        # all_peptide_forest_percolator.append(filtered_peptide_forest_perc)


        # # Peptide Forest for Open Mod

        # uc.params['prefix'] = 'open_pf_only_' + mzml_basename.replace('.mzML', '')
        # uc.params['peptide_forest_initial_engine'] = 'msfragger_2_3'
        # open_validated_peptide_forest_alone = uc.validate(
        #     input_file=open_unvalidated,
        #     engine='peptide_forest',
        # )
        # uc.params['csv_filter_rules'] = [
        #     ['Is decoy', 'equals', 'false'],
        #     ['q-value_RF-reg','lte', 0.01],
        #     ['Conflicting uparam', 'contains_not', 'enzyme'],
        # ]
        # open_filtered_peptide_forest_alone = uc.execute_misc_engine(
        #     input_file = open_validated_peptide_forest_alone,
        #     engine='filter_csv',
        # )
        # all_open_peptide_forest_alone.append(open_filtered_peptide_forest_alone)

        # uc.params['prefix'] = 'open_pf_perc_' + mzml_basename.replace('.mzML', '')
        # open_validated_peptide_forest_perc = uc.validate(
        #     input_file=open_validated,
        #     engine='peptide_forest',
        # )
        # uc.params['csv_filter_rules'] = [
        #     ['Is decoy', 'equals', 'false'],
        #     ['q-value_RF-reg','lte', 0.01],
        #     ['Conflicting uparam', 'contains_not', 'enzyme'],
        # ]
        # open_filtered_peptide_forest_perc = uc.execute_misc_engine(
        #     input_file = open_validated_peptide_forest_perc,
        #     engine='filter_csv',
        # )
        # all_open_peptide_forest_percolator.append(open_filtered_peptide_forest_perc)
        # uc.params['prefix'] = ''

        # # validated_peptide_forest_de_novo = uc.validate(
        # #     input_file=,
        # #     engine='peptide_forest',
        # # )
        # # all_peptide_forest_de_novo.append(validated_peptide_forest_de_novo)

    ##################################
    # Merge section
    ##################################

    # uc.params['prefix'] = 'closed_search_perc_'
    # merged_all_perc = uc.execute_misc_engine(
    #     input_file = all_percolator_results,
    #     engine='merge_csvs',
    #     merge_duplicates=True,
    # )

    # uc.params['prefix'] = 'closed_search_combined_'
    # merged_all_combined_pep = uc.execute_misc_engine(
    #     input_file = all_combined_pep_results,
    #     engine='merge_csvs',
    #     merge_duplicates=True,
    # )

    uc.params['prefix'] = 'closed_search_'
    merged_all_peptide_forest_alone = uc.execute_misc_engine(
        input_file = all_peptide_forest_alone,
        engine='merge_csvs',
        merge_duplicates=True,
    )

    # uc.params['prefix'] = 'closed_search_pf_perc_'
    # merged_all_peptide_forest_perc = uc.execute_misc_engine(
    #     input_file = all_peptide_forest_percolator,
    #     engine='merge_csvs',
    #     merge_duplicates=True,
    # )

    # # merged_all_peptide_forest_de_novo = uc.execute_misc_engine(
    # #     input_file = all_peptide_forest_de_novo,
    # #     engine='merge_csvs',
    # #     merge_duplicates=True,
    # # )

    # uc.params.update({
    #     'visualization_column_names': ['Modifications', 'Sequence'],
    #     'visualization_label_positions':{
    #         '0': 'Percolator',
    #         '1': 'Combined PEP',
    #         '2': 'Peptide Forest',
    #         '3': 'Peptide Forest + Percolator'
    #     }
    # })
    # uc.visualize(
    #     input_files=[
    #         merged_all_perc,
    #         merged_all_combined_pep,
    #         merged_all_peptide_forest_alone,
    #         merged_all_peptide_forest_perc, 
    #     ],
    #     engine='venndiagram_1_1_0',
    #     output_file_name='closed_search_combined_venndiagram_peptides_filtered.svg'
    # )

    # uc.params['visualization_column_names'] = ['Spectrum Title', 'Sequence', 'Modifications']
    # uc.visualize(
    #     input_files=[
    #         merged_all_perc,
    #         merged_all_combined_pep,
    #         merged_all_peptide_forest_alone,
    #         merged_all_peptide_forest_perc, 
    #     ],
    #     engine='venndiagram_1_1_0',
    #     output_file_name='closed_search_combined_venndiagram_psms_filtered.svg'
    # )

    # # Open Mod validation

    # open_merged_all_perc = uc.execute_misc_engine(
    #     input_file = all_open_validated,
    #     engine='merge_csvs',
    #     merge_duplicates=True,
    # )

    # open_merged_all_combined_pep = uc.execute_misc_engine(
    #     input_file = all_open_combined_pep,
    #     engine='merge_csvs',
    #     merge_duplicates=True,
    # )

    # open_merged_all_peptide_forest_only = uc.execute_misc_engine(
    #     input_file = all_open_peptide_forest_alone,
    #     engine='merge_csvs',
    #     merge_duplicates=True,
    # )

    # open_merged_all_peptide_forest_perc = uc.execute_misc_engine(
    #     input_file = all_open_peptide_forest_percolator,
    #     engine='merge_csvs',
    #     merge_duplicates=True,
    # )

    # uc.params.update({
    #     'visualization_column_names': ['Modifications', 'Sequence'],
    #     'visualization_label_positions':{
    #         '0': 'Percolator',
    #         '1': 'Combined PEP',
    #         '2': 'Peptide Forest',
    #         '3': 'Peptide Forest + Percolator'
    #     }
    # })
    # uc.visualize(
    #     input_files=[
    #         open_merged_all_perc,
    #         open_merged_all_combined_pep,
    #         open_merged_all_peptide_forest_only,
    #         open_merged_all_peptide_forest_perc, 
    #     ],
    #     engine='venndiagram_1_1_0',
    #     output_file_name='open_search_combined_venndiagram_peptides.svg'
    # )

    # uc.params['visualization_column_names'] = ['Spectrum Title', 'Sequence', 'Modifications']
    # uc.visualize(
    #     input_files=[
    #         open_merged_all_perc,
    #         open_merged_all_combined_pep,
    #         open_merged_all_peptide_forest_only,
    #         open_merged_all_peptide_forest_perc, 
    #     ],
    #     engine='venndiagram_1_1_0',
    #     output_file_name='open_search_combined_venndiagram_psms.svg'
    # )

    # # uc.params.update({
    # #     'mgf_input_files_list' : all_mgf_files,
    # #     'validation_score_field': 'PEP',
    # #     'bigger_scores_better': False,
    # # })

    # # ptminer_all_perc = uc.validate(
    # #     input_file=open_merged_all_perc,
    # #     engine='ptminer',
    # #     # force=True,
    # # )

    # # uc.params["csv_filter_rules"] = [
    # #     ["Is decoy", "equals", "false"],
    # #     ["PEP", "lte", 0.01],
    # #     ["Conflicting uparam", "contains_not", "enzyme"],
    # # ]

    # # filtered_ptminer_all_perc = uc.execute_misc_engine(
    # #     input_file=ptminer_all_perc,
    # #     engine='filter_csv_1_0_0',
    # # )

    # # ptminer_all_combined_pep = uc.validate(
    # #     input_file=open_merged_all_combined_pep,
    # #     engine='ptminer',
    # #     # force=True,
    # # )

    # # uc.params['csv_filter_rules'] = [
    # #     ['Is decoy', 'equals', 'false'],
    # #     ['combined PEP','lte', 0.01],
    # #     ['Conflicting uparam', 'contains_not', 'enzyme'],
    # # ]

    # # filtered_ptminer_all_combined_pep = uc.execute_misc_engine(
    # #     input_file=ptminer_all_combined_pep,
    # #     engine='filter_csv_1_0_0',
    # # )








    # # params = {optimized_from_sweep, wide_search_arams}
    # # for file in files:
    # #     percolator_results = []
    # #     unvalidated_results = []
    # #     for open_mod_engine in [
    # #         'msfragger',
    # #         'moda',
    # #         'pipi',
    # #         'taggraph',
    # #     ]:
    # #         unified_search_results = uc.search(file, prot_db_engine)
    # #         unvalidated_results.append(unified_search_results)

    # #         vd_engine = 'percolator_3_4'
    # #         validated_results = uc.valdiate(vd_engine, unified_search_results)
    # #         percolator_results.append(validated_results)

    # #     vd_engine = 'peptide_forest'
    # #     peptide_forest_alone = uc.valdiate(vd_engine, unvalidated_results)
    # #     peptide_forest_percolator = uc.valdiate(vd_engine, percolator_results)

    # #     venn_diagram = uc.visualize([percolator_results, peptide_forest_alone, peptide_forest_percolator])


if __name__ == '__main__':
    main(
        folder=sys.argv[1],
        database=sys.argv[2],
    )


