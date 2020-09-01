import os
import json

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONF_DIR = os.path.join(BASE_DIR, 'conf')
DATA_DIR = os.path.join(BASE_DIR, 'data')

HPO_VERSION = '2020-03-13'

HPO_DIR = os.path.join(DATA_DIR, 'hpo', HPO_VERSION)
HPO_OBO_PATH = os.path.join(HPO_DIR, 'hp.obo')
HPO_PHENOTYPE_TO_GENE_PATH = os.path.join(HPO_DIR, 'phenotype_to_genes.txt')
HPO_GENE_TO_PHENOTYPE_PATH = os.path.join(HPO_DIR, 'genes_to_phenotype.txt')
HPO_GENE_SIM_PATH = os.path.join(HPO_DIR, 'hpo-gene-similarity.pkl')
HPO_PARAMS_JSON_FILE = os.path.join(CONF_DIR, 'psea_params_full.json')
HPO_PARAMS = json.load(open(HPO_PARAMS_JSON_FILE))

try:
    from psea.local_settings import *
except ImportError:
    pass