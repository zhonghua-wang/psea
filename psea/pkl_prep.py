import itertools
import itertools as it
import pickle

import pandas as pd

from hpo import HPO
from tqdm import tqdm


def create_hpo_hpo_pkl(hpo_config: dict, use_full: bool, seri_full: bool):
    hpo = HPO()
    result_dict = {}
    for hpo1, hpo2 in tqdm(
            itertools.combinations_with_replacement(hpo.phen2gene_hpo_list, 2),
            desc='building HPO-HPO similarity file'
    ):
        similarity_func = hpo.get_eric_similarity_unweighted
        result_dict[(hpo1, hpo2)] = similarity_func(hpo1, hpo2)
    with open(hpo_config['hpo_hpo_pkl'], 'wb') as f:
        pickle.dump(result_dict, f)


def create_hpo_gene_pkl(hpo_config: dict, use_full: bool, seri_full: bool, weighted=False):
    hpo = HPO()

    def _func(hpo_id):
        _dict = {}
        for gene_id, _hpo_li in hpo.gene_hpo_series.items():
            _dict[gene_id] = hpo.get_hpo_set_eric_similarity(
                [hpo_id], _hpo_li)
        _seri = pd.Series(_dict)
        _seri.name = hpo_id
        return _seri

    def _func_weighted(hpo_id):
        _dict = {}
        for gene_id, _df in hpo.gene_phen_map_df.sort_values(['gene_id', 'freq_weight'], ascending=False).groupby(
                'gene_id'):
            _df = _df.drop_duplicates(['hpo_id'], keep='first')
            _dict[gene_id] = (
                    _df['hpo_id'].apply(lambda x: hpo.get_eric_similarity_unweighted(
                        hpo_id, x)) * _df['freq_weight']
            ).sum()
        _seri = pd.Series(_dict)
        _seri.name = hpo_id
        return _seri

    func = _func_weighted if weighted else _func

    _df = pd.Series(hpo.phen2gene_hpo_list,
                    hpo.phen2gene_hpo_list).parallel_apply(func)
    # return _df
    with open(hpo_config['gene_hpo_pkl'], 'wb') as w:
        pickle.dump(_df, w)
