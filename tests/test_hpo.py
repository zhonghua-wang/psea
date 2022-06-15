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
