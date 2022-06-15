import pandas as pd
import obonet
import networkx as nx
import math
import os

from psea.hpo_config import HPOConfig
import collections as co
import functools as ft
import itertools as it


class HPO:
    def __init__(
            self, use_full=True
            # self, use_full=True, version=settings.HPO_VERSION,
            # obo_path=settings.HPO_OBO_PATH, gene2phen_path=settings.HPO_GENE_TO_PHENOTYPE_PATH,
            # phen2gene_path=settings.HPO_PHENOTYPE_TO_GENE_PATH, gene_hpo_pkl=settings.HPO_GENE_SIM_PATH,
            # hpo_hpo_pkl=None, old=False,
    ):
        self.hpo_conf = HPOConfig()
        self.version = self.hpo_conf.version
        self._ancestor_dict = {}
        self.full_graph = obonet.read_obo(self.hpo_conf.obo_path)
        self.hpo_set = nx.ancestors(self.full_graph, 'HP:0000118')
        self.hpo_set.add('HP:0000118')
        self.graph = self.full_graph.subgraph(self.hpo_set)

        self.gene2phen = pd.read_csv(
            self.hpo_conf.genes_to_phenotype_path, skiprows=[0],
            names=['gene_id', 'gene_name', 'hpo_id', 'hpo_name', 'freq_raw', 'freq', 'info', 'source', 'source_id'],
            sep='\t'
        )

        self.phen2gene = pd.read_csv(
            self.hpo_conf.phenotype_to_genes_path, skiprows=[0],
            names=['hpo_id', 'hpo_name', 'gene_id', 'gene_name', 'info', 'source', 'source_id'], sep='\t'
        )
        self.phen2gene = self.phen2gene.loc[self.phen2gene['hpo_id'].isin(
            self.hpo_set)]
        self.gene2phen = self.gene2phen.loc[self.gene2phen['hpo_id'].isin(
            self.hpo_set)]
        if use_full:
            self.gene_phen_map_df = self.phen2gene
        else:
            self.gene_phen_map_df = self.gene2phen
        self.gene_set = set(self.phen2gene['gene_id'])
        self.gene_ic_dict = dict(
            [
                (k, -math.log(v / len(self.gene_phen_map_df)))
                for k, v in co.Counter(self.gene_phen_map_df['gene_id']).most_common()
            ]
        )
        self.hpo_ic_dict = dict(
            [
                (x[0], -math.log(x[1] / len(self.gene_phen_map_df)))
                for x in co.Counter(self.gene_phen_map_df['hpo_id']).most_common()
            ]
        )
        self.eric_ic_dict = dict(
            (hpo_id, -math.log(len(gene_list) / len(self.gene_set)))
            for hpo_id, gene_list in self.gene_phen_map_df.groupby('hpo_id')
        )

        self.phen2gene_hpo_list = list(set(self.phen2gene['hpo_id']))
        self.gene_hpo_series = self.gene_phen_map_df.groupby('gene_id')[
            'hpo_id'].apply(list)
        self.gene_hpo_size = self.gene_hpo_series.map(len)
        self.gene_hpo_size.name = 'size'

        self.gene_hpo_sim_df = None
        self.hpo_hpo_sim_df = None

        # if gene_hpo_pkl is None or not os.path.exists(gene_hpo_pkl):
        #     self.gene_hpo_sim_df = None
        # else:
        #     with open(gene_hpo_pkl, 'rb') as f:
        #         self.gene_hpo_sim_df = pickle.load(f)
        # if hpo_hpo_pkl is None or not os.path.exists(hpo_hpo_pkl):
        #     self.hpo_hpo_sim_pkl = None
        # else:
        #     with open(hpo_hpo_pkl, 'rb') as f:
        #         self.hpo_hpo_sim_pkl = pickle.load(f)

        self.gene_df = self.gene2phen.groupby('gene_id')[['gene_name']].first()

    # def get_frequency(self, hpo_id, gene_id, include_child=True):
    #     if self._freq_dict is None:
    #         self._freq_dict = self.gene2phen.groupby(['hpo_id', 'gene_id'])['freq_weight'].agg(max).to_dict()
    #     weight = self._freq_dict.get((hpo_id, gene_id))
    #     if weight is not None or not include_child:
    #         return weight
    #
    #     weight_li = []
    #     for child in self.get_descendant_node_id_set(hpo_id):
    #         _weight = self._freq_dict.get((child, gene_id), 0)
    #         if _weight == HPOFrequency.max_weight():
    #             return _weight
    #         weight_li.append(_weight)
    #     return max(weight_li) or HPOFrequency.DEFAULT_WEIGHT

    def get_ancestor_node_id_set(self, hpo_id):
        """get ancestors hpo id set of specified HPO ID

        Arguments:
            hpo_id {str} -- HPO term id, eg. HP:0000118
        """
        return nx.descendants(self.graph, hpo_id)

    @property
    def ancestor_dict(self):
        if self._ancestor_dict == {}:
            for hpo in self.phen2gene_hpo_list:
                self._ancestor_dict[hpo] = self.get_ancestor_node_id_set(hpo)
        return self._ancestor_dict

    def get_ancestor_hpo_set(self, hpo_id):
        """
        get ancestor hpo id set of specified HPO ID
        :param hpo_id: HPO ID, i.e. HP:0000525
        :return: HPO ID set
        """
        comm = self.ancestor_dict.get(hpo_id)
        if comm is None:
            return {"HP:0000118"}
        comm.add("HP:0000118")
        comm.add(hpo_id)
        return comm

    def get_node(self, hpo_id: str, full=False):
        """get hpo entry by HPO ID

        Arguments:
            hpo_id {str} -- HPO term id, eg. HP:0000118
        """

        if full:
            return self.full_graph.nodes[hpo_id]
        return self.graph.nodes[hpo_id]

    def get_descendant_node_id_set(self, hpo_id):
        """get descendants node id set of specified HPO ID

        Arguments:
            hpo_id {str} -- HPO term id, eg. HP:0000118
        """
        return nx.ancestors(self.graph, hpo_id)

    def get_common_ancestors(self, *hpo_id_list):
        # TODO: get latest ancestor
        return ft.reduce(
            lambda x, y: x.intersection(y),
            [self.get_ancestor_hpo_set(hpo) for hpo in hpo_id_list]
        )

    def get_ic(self, hpo_id):
        return self.hpo_ic_dict.get(hpo_id)

    def get_eric_ic(self, hpo_id):
        return self.eric_ic_dict.get(hpo_id, 0)

    def get_ic_similarity(self, hpo_1, hpo_2):
        sim_list = [self.get_ic(x)
                    for x in self.get_common_ancestors(hpo_1, hpo_2)]
        sim_list = list(filter(lambda x: x is not None, sim_list))
        return max(sim_list) if len(sim_list) > 0 else 0

    def get_hpo_sim_pkl(self, hpo1, hpo2):
        return self.hpo_hpo_sim_pkl.get((hpo1, hpo2), 0) or self.hpo_hpo_sim_pkl.get((hpo2, hpo1), 0)

    def get_hpo_set_sim_pkl(self, hpo_set_1, hpo_set_2, cutoff=0):
        return sum(
            filter(lambda x: x >= cutoff,
                   [self.get_hpo_sim_pkl(x, y) for x, y in it.product(hpo_set_1, hpo_set_2)])
        )

    def get_set_ic_sim(self, set1, set2):
        """Calculate IC similarity (rawscore) of the specified two HPO sets

        Arguments:
            set1 {set} -- first HPO set
            set2 {set} -- second HPO set

        Returns:
            float -- IC similarity value
        """

        return sum([self.get_ic_similarity(x1, x2) for x1, x2 in it.product(set1, set2)])

    def get_parent(self, hpo_id, valid_only=False):
        """
        :param hpo_id: HPO ID, i.e. HP:0000365
        :param valid_only:
        :return:
        """
        if hpo_id == 'HP:0000001':
            return []
        if valid_only:
            try:
                parent_list = self.get_node(hpo_id).get('is_a')
            except KeyError:
                return []
            valid_li = []
            for parent_id in parent_list:
                if parent_id in self.eric_ic_dict.keys():
                    valid_li.append(parent_id)
            if len(valid_li) > 0:
                return valid_li
            else:
                for parent_id in parent_list:
                    valid_li.extend(self.get_parent(
                        parent_id, valid_only=valid_only))
            return valid_li
        try:
            return self.get_node(hpo_id).get('is_a')
        except KeyError:
            return []

    def get_parent_set(self, hpo_id):
        return set(self.graph.successors(hpo_id))

    def get_sibling_set(self, hpo_id):
        sibling_set = set(
            list(
                it.chain.from_iterable(
                    [self.get_children_set(x) for x in self.get_parent(hpo_id)]
                )
            )
        )
        sibling_set.remove(hpo_id)
        return sibling_set

    def get_children_set(self, hpo_id):
        return set(self.graph.predecessors(hpo_id))

    def get_eric_similarity_unweighted(self, hpo_id_1, hpo_id_2):
        return max(
            0,
            2 * max([self.get_eric_ic(x) for x in self.get_common_ancestors(hpo_id_1, hpo_id_2)]) - min(
                self.get_eric_ic(hpo_id_1), self.get_eric_ic(hpo_id_2))
        )

    def get_hpo_gene_similarity(self, hpo_id, gene_id):
        try:
            return self.gene_hpo_sim_df.loc[hpo_id][gene_id]
        except KeyError:
            return 0

    def get_hpo_set_eric_similarity(self, hpo_set_1, hpo_set_2, cutoff=0):
        return sum(
            filter(lambda x: x >= cutoff,
                   [self.get_eric_similarity_unweighted(x, y) for x, y in it.product(hpo_set_1, hpo_set_2)])
        )

    def gene_symbol_to_hpo(self, gene_symbol, full=False):
        if full:
            tmp_df = self.phen2gene[self.phen2gene['gene_name'] == gene_symbol]
        else:
            tmp_df = self.gene2phen[self.gene2phen['gene_name'] == gene_symbol]
        return tmp_df['hpo_id'].to_list()

    def hpo_to_gene_symbol(self, hpo_id, full=False):
        tmp_df = self.phen2gene[self.phen2gene['hpo_id'] == hpo_id] if full \
            else self.gene2phen[self.gene2phen['hpo_id'] == hpo_id]

        return tmp_df['gene_name'].to_list()

    def valid_hpo(self) -> pd.DataFrame:
        hpo_tmp = self.graph.copy()
        invalid_hpo = set(self.graph.nodes).difference(
            set(self.phen2gene_hpo_list))
        for hpo_el in invalid_hpo:
            hpo_tmp.remove_node(hpo_el)
        result_df = pd.DataFrame()
        for hpo_el in hpo_tmp.nodes:
            node = hpo_tmp.nodes[hpo_el]
            node['hpo_id'] = hpo_el
            result_df = result_df.append(node, ignore_index=True)
        return result_df

    def get_parents_id(self, hpo_id, full=True):
        hpo = self.get_node(hpo_id, full)
        if hpo.get('is_a') is None:
            return []
        li = []
        for parent_id_el in hpo['is_a']:
            li.append(parent_id_el)
            li.extend(self.get_parents_id(parent_id_el, full=True))
            # li.reverse()
        return li

    def ordered_hpo_list(self, full=False):
        li = []
        if full:
            graph = self.full_graph
        else:
            graph = self.graph
        for hpo_id in graph.nodes:
            li.append(hpo_id)
            li.extend(self.get_parents_id(hpo_id, full))
        li.reverse()
        return list(co.OrderedDict.fromkeys(li).keys())

    def get_valid_hpo(self, hpo_id):
        return self.get_parent(hpo_id, valid_only=True)

    def get_hpo_resnik_similarity(self, hpo1: str, hpo2: str):
        return max([self.get_eric_ic(x) for x in self.get_common_ancestors(hpo1, hpo2)])

    def get_hpo_set_resnik_similarity(self, hpo_set1: list, hpo_set2: list, cutoff=0):
        return sum(
            filter(lambda x: x >= cutoff,
                   [self.get_hpo_resnik_similarity(x, y) for x, y in it.product(hpo_set1, hpo_set2)])
        )
