import os
import json

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONF_DIR = os.path.join(BASE_DIR, 'conf')
# DATA_DIR = os.path.join(BASE_DIR, 'data')
LOCAL_DATA_DIR = os.path.join(BASE_DIR, '.data')
HPO_DIR = os.path.join(LOCAL_DATA_DIR, 'hpo')
HPO_LOCAL_CONF_FILE = os.path.join(LOCAL_DATA_DIR, 'hpo.json')

# HPO_PHENOTYPE_TO_GENE_PATH = os.path.join(HPO_DIR, 'phenotype_to_genes.txt')
# HPO_GENE_TO_PHENOTYPE_PATH = os.path.join(HPO_DIR, 'genes_to_phenotype.txt')
# HPO_GENE_SIM_PATH = os.path.join(HPO_DIR, 'hpo-gene-similarity.pkl')
HPO_PARAMS_JSON_FILE = os.path.join(CONF_DIR, 'psea_params_full.json')
HPO_PARAMS = json.load(open(HPO_PARAMS_JSON_FILE))
#
# GENES_TO_PHENOTYPE_URL = 'http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt'
# PHENOTYPE_TO_GENE_URL = 'http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt'
#
# PHENOTYPE_ANNOTATION_URL = 'http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa'
# PHENOTYPE_ONTOLOGY_URL = 'http://purl.obolibrary.org/obo/hp.obo'

try:
    from psea.local_settings import *
except ImportError:
    pass

