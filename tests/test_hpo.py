import os
from psea.hpo import HPO
from psea.hpo_config import HPOConfig


def test_hpo_config():
    hpo_conf = HPOConfig()
    assert os.path.exists(hpo_conf.obo_path)
    assert os.path.exists(hpo_conf.phenotype_to_genes_path)
    assert os.path.exists(hpo_conf.genes_to_phenotype_path)
    assert os.path.exists(hpo_conf.annotation_path)
    assert hpo_conf.version.startswith('20')


def test_hpo():
    hpo = HPO()
    assert hpo.version.startswith('20')
    assert len(hpo.get_parent(hpo_id='HP:0000001')) == 0
    assert len(hpo.gene_hpo_sim_df) >= 10000
    assert hpo.get_hpo_sim_pkl('HP:0000360', 'HP:0000365') > 0


def test_psea():
    from psea.PSEA import PSEA
    psea = PSEA()
    _df = psea.predict_pkl(['HP:0000360', 'HP:0000365'])
    _dfv = psea.predict_pkl_verbose(['HP:0000360', 'HP:0000365'])
    assert psea.hpo.version.startswith('20')
    assert len(_df) >= 10
    assert 9132 in _df.index
    print(_df.head())
    print(_dfv.head())
