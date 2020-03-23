params = {optimized_from_sweep}
uc = ursgal.UController(params=params, profile='QExactive+', verbose=False)

for file in files:
    percolator_results = []
    unvalidated_results = []
    for prot_db_engine in [
        'xtandem',
        'msgfplus',
        'msfragger',
        'mascot',
        'omssa',
        'msamanda',
    ]:
        unified_search_results = uc.search(file, prot_db_engine)
        unvalidated_results.append(unified_search_results)

        vd_engine = 'percolator_3_4'
        validated_results = uc.valdiate(vd_engine, unified_search_results)
        percolator_results.append(validated_results)

    de_novo_results = []
    for de_novo_engine in [
        pnovo,
        novor,
        deepnovo
    ]:
        unified_search_results = uc.search(file, prot_db_engine)
        de_novo_results.append(unified_search_results)

    vd_engine = 'peptide_forest'
    peptide_forest_alone = uc.valdiate(vd_engine, unvalidated_results)
    peptide_forest_alone = uc.valdiate(vd_engine, de_novo_results)
    peptide_forest_percolator = uc.valdiate(vd_engine, percolator_results)

    venn_diagram = uc.visualize([percolator_results, peptide_forest_alone, peptide_forest_percolator])

params = {optimized_from_sweep, wide_search_arams}
for file in files:
    percolator_results = []
    unvalidated_results = []
    for open_mod_engine in [
        'msfragger',
        'moda',
        'pipi',
        'taggraph',
    ]:
        unified_search_results = uc.search(file, prot_db_engine)
        unvalidated_results.append(unified_search_results)

        vd_engine = 'percolator_3_4'
        validated_results = uc.valdiate(vd_engine, unified_search_results)
        percolator_results.append(validated_results)

    vd_engine = 'peptide_forest'
    peptide_forest_alone = uc.valdiate(vd_engine, unvalidated_results)
    peptide_forest_percolator = uc.valdiate(vd_engine, percolator_results)

    venn_diagram = uc.visualize([percolator_results, peptide_forest_alone, peptide_forest_percolator])


