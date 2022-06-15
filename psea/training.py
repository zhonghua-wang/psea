import itertools
import pickle
import pandas as pd
from psea.hpo import HPO
from tqdm import tqdm
from loguru import logger

tqdm.pandas()


def create_hpo_hpo_pkl(hpo):
    result_dict = {}
    for hpo1, hpo2 in tqdm(itertools.combinations_with_replacement(
            hpo.phen2gene_hpo_list, 2), total=len(hpo.phen2gene_hpo_list) * (len(hpo.phen2gene_hpo_list) + 1) / 2
    ):
        similarity_func = hpo.get_eric_similarity_unweighted
        result_dict[(hpo1, hpo2)] = similarity_func(hpo1, hpo2)
    with open(hpo.hpo_conf.hpo_hpo_similarity_pkl_path, 'wb') as f:
        pickle.dump(result_dict, f)


def create_hpo_gene_pkl(hpo: HPO, nb_workers=1):
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

    func = _func

    if nb_workers > 1:
        try:
            from pandarallel import pandarallel
            pandarallel.initialize(nb_workers=nb_workers)
            _df = pd.Series(hpo.phen2gene_hpo_list,
                            hpo.phen2gene_hpo_list).parallel_apply(func)
        except BaseException:
            Warning('Failed to run in parallel mode, using 1 worker')
            __df = pd.Series(hpo.phen2gene_hpo_list, hpo.phen2gene_hpo_list).progress_apply(func)
    else:
        _df = pd.Series(hpo.phen2gene_hpo_list, hpo.phen2gene_hpo_list).progress_apply(func)
    # return _df
    with open(hpo.hpo_conf.hpo_gene_similarity_pkl_name, 'wb') as w:
        pickle.dump(_df, w)


def build(nb_workers=1):
    hpo = HPO()
    logger.info(f'Building HPO-HPO file for HPO.{hpo.version}')
    create_hpo_hpo_pkl(hpo)
    logger.info(f'Building Gene-HPo file')
    create_hpo_gene_pkl(hpo, nb_workers=nb_workers)
    logger.info('Finish building HPO data file')

# if __name__ == '__main__':
#     pandarallel.initialize(progress_bar=True, nb_workers=4)
#     for use_full, seri_full in it.product([True, False], repeat=2):
#         version = f'2021-12-17-{use_full}-{seri_full}'
#         hpo_config = hpo_util.get_hpo_config(version)
#         logger.info(f'Creating... HPO-HPO pickle file for version {version}')
#         # create_hpo_hpo_pkl(hpo_config, use_full, seri_full=seri_full)
#         logger.info(f'finished creating HPO-HPO pickle file')
#         create_hpo_gene_pkl(hpo_config, use_full=use_full, seri_full=seri_full)
#         logger.info(f'finished creating HPO-Gene pickle file')
