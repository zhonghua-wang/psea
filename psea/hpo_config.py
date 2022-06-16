import json
from datetime import datetime
import os
from urllib import request
from tqdm import tqdm
from psea import settings
from dataclasses import dataclass


@dataclass
class HPOConfig:
    """
    HPO config
    """

    _obo_name: str = 'hp.obo'
    _genes_to_phenotype_name: str = 'genes_to_phenotype.txt'
    _phenotype_to_genes_name: str = 'phenotype_to_genes.txt'
    _annotation_name: str = 'phenotype.hpoa'
    _hpo_hpo_similarity_pkl_name: str = 'hpo-hpo-similarity.pkl'
    _hpo_gene_similarity_pkl_name: str = 'hpo-gene-similarity.pkl'
    _version: str = ''

    GENES_TO_PHENOTYPE_URL = 'http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt'
    PHENOTYPE_TO_GENE_URL = 'http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt'

    PHENOTYPE_ANNOTATION_URL = 'http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa'
    PHENOTYPE_ONTOLOGY_URL = 'http://purl.obolibrary.org/obo/hp.obo'

    @staticmethod
    def my_hook(t):
        last_b = [0]

        def inner(b=1, bsize=1, tsize=None):
            if tsize is not None:
                t.total = tsize
            if b > 0:
                t.update((b - last_b[0]) * bsize)
            last_b[0] = b

        return inner

    def init_hpo(self, override=False):
        self._version = datetime.now().strftime("%Y-%m-%d")
        # if os.path.exists(self.hpo_dir) and not override:
        #     print(
        #         f'HPO data directory ({self.hpo_dir}) already exists. Skipping download. Use --override to force '
        #         f'download.')
        #     return
        os.makedirs(self.hpo_dir, exist_ok=True)
        print('1. Downloading HPO ontology data...')
        if os.path.exists(self.obo_path) and not override:
            print(
                f'HPO data ({self.hpo_dir}) already exists. Skipping download. Use --override to force '
                f'download.')
            return
        with tqdm(unit='B', unit_scale=True, miniters=1) as t:
            # filename, _ = urlretrieve(url, filename=destination, reporthook=my_hook(t))
            request.urlretrieve(self.PHENOTYPE_ONTOLOGY_URL, self.obo_path, reporthook=HPOConfig.my_hook(t))
        print('2. Downloading HPO annotation data...')
        with tqdm(unit='B', unit_scale=True, miniters=1) as t:
            request.urlretrieve(self.PHENOTYPE_ANNOTATION_URL, self.annotation_path, reporthook=HPOConfig.my_hook(t))
        print('3. Downloading HPO genes to phenotype data...')
        with tqdm(unit='B', unit_scale=True, miniters=1) as t:
            request.urlretrieve(self.GENES_TO_PHENOTYPE_URL, self.genes_to_phenotype_path,
                                reporthook=HPOConfig.my_hook(t))
        print('4. Downloading HPO phenotype to genes data...')
        with tqdm(unit='B', unit_scale=True, miniters=1) as t:
            request.urlretrieve(self.PHENOTYPE_TO_GENE_URL, self.phenotype_to_genes_path,
                                reporthook=HPOConfig.my_hook(t))

        with open(settings.HPO_LOCAL_CONF_FILE, 'w') as f:
            json.dump({'version': self._version}, f)
        print(f'HPO version was updated to date: {self._version}')

    @property
    def version(self):
        if self._version == '':
            if not os.path.exists(settings.HPO_LOCAL_CONF_FILE):
                raise FileNotFoundError(f'HPO data file not found, please download it first')
            conf = json.load(open(settings.HPO_LOCAL_CONF_FILE))
            self._version = conf['version']
        return self._version

    @property
    def hpo_dir(self):
        # return os.path.join(settings.HPO_DIR, self.version)
        return settings.HPO_DIR

    @property
    def obo_path(self):
        return os.path.join(self.hpo_dir, self._obo_name)

    @property
    def annotation_path(self):
        return os.path.join(self.hpo_dir, self._annotation_name)

    @property
    def genes_to_phenotype_path(self):
        return os.path.join(self.hpo_dir, self._genes_to_phenotype_name)

    @property
    def phenotype_to_genes_path(self):
        return os.path.join(self.hpo_dir, self._phenotype_to_genes_name)

    @property
    def hpo_hpo_similarity_pkl_path(self):
        # _path = os.path.join(self.hpo_dir, self._hpo_hpo_similarity_pkl_name)
        # if not os.path.exists(_path):
        #     raise FileNotFoundError(
        #         f'HPO similarity file not found, please run HPOConfig.compute_hpo_hpo_similarity() first')
        return os.path.join(self.hpo_dir, self._hpo_hpo_similarity_pkl_name)

    @property
    def hpo_gene_similarity_pkl_path(self):
        # _path = os.path.join(self.hpo_dir, self._hpo_gene_similarity_pkl_name)
        # if not os.path.exists(_path):
        #     raise FileNotFoundError(
        #         f'HPO gene similarity file not found, please run HPOConfig.compute_hpo_gene_similarity() first')
        return os.path.join(self.hpo_dir, self._hpo_gene_similarity_pkl_name)
