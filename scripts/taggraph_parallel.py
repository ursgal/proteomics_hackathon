#!/usr/bin/env python3.4
# encoding: utf-8

import ursgal
import os
import sys
import shutil
import glob
import multiprocessing


def run_taggraph(ms_file):
    uc = ursgal.UController(
        profile='QExactive+',
        params={
            'database': '../../data/Uniprot_swissprot_TREMBL_cRAP_target_decoy.fasta',
            'modifications': [
                'M,opt,any,Oxidation',        # Met oxidation
                'C,fix,any,Carbamidomethyl',  # Carbamidomethylation
                'N,opt,any,Deamidated',
                'Q,opt,any,Deamidated',
            ],
            'peptide_mapper_class_version': 'UPeptideMapper_v3',
            '-xmx': '4500m',
            'cpus': 11,
            'csv_filter_rules' : [
                ['Is decoy', 'equals', 'false'],
                ['PEP','lte', 0.01],
                ['Conflicting uparam', 'contains_not', 'enzyme'],
            ],
            'enzyme': 'trypsin',
            'pymzml_spec_id_attribute': {'ID': None},
            'frag_mass_tolerance'       : 20,
            'frag_mass_tolerance_unit'  : 'ppm',
            'precursor_mass_tolerance_plus' : 5,
            'precursor_mass_tolerance_minus' : 5,
        }
    )

    engine = 'deepnovo_v2'
    unified_file_list = []
    'de_novo_TN_CSF_062617_02___unified_merged.csv'
    dirname = '/media/external/Projects/proteomics_blog_hackathon/data/Q-Exactive_Plus/merged_denovo/de_novo_results'
    # dirname = f'/media/external/Projects/proteomics_blog_hackathon/data/Q-Exactive_Plus/merged_denovo/de_novo_results/{engine}'
    ms_base = os.path.basename(ms_file.replace('.mzML', ''))
    unified_search_result = os.path.abspath(os.path.join(
        dirname,
        f'de_novo_{ms_base}___unified_merged.csv',
        # f'{ms_base}_{engine}_unified.csv'
    ))
    assert os.path.isfile(unified_search_result)
    merged_files = unified_search_result
    uc.params.update({'modifications': [
                'M,opt,any,Oxidation',        # Met oxidation
                'C,fix,any,Carbamidomethyl',  # Carbamidomethylation
                'N,opt,any,Deamidated',
                'Q,opt,any,Deamidated',
            ],})
    uc.params['de_novo_results'] = merged_files
    taggraph_search_result = uc.search(
        input_file=ms_file,
        engine='tag_graph_1_8_0',
    )
    return


def main(database=None, ms_folder=None):
    '''
    Executes a search with OMSSA, XTandem and MS-GF+ on the BSA1.mzML
    input_file

    usage:
        ./simple_example_search.py

    Note:
        Myrimatch does not work with this file.
        To use MSAmanda on unix platforms, please install mono
        (http://www.mono-project.com/download)

    '''

    ms_files = [f for f in glob.glob(os.path.join(ms_folder, "*.mzML"))]

    pool = multiprocessing.Pool(2)
    pool.map(run_taggraph, ms_files)
    return


if __name__ == '__main__':
    main(ms_folder=sys.argv[1], database=sys.argv[2])
