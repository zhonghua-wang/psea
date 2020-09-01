import pandas as pd
import itertools as it
from psea import settings
from psea.hpo import HPO
import warnings


class PSEA(object):
    def __init__(self, hpo: HPO = None):
        if hpo is None:
            self.hpo = HPO()
        else:
            self.hpo = hpo
        self.param = settings.HPO_PARAMS[self.hpo.version]

    def size_rs_mean(self, size):
        return self.param['MEAN_PARAM'] * size

    def size_rs_std(self, size):
        return self.param['STD_PARAM_1'] * (size ** self.param['STD_PARAM_2'])

    def _zscore(self, rs, size):
        if size == 0:
            return -999
        return (rs - self.size_rs_mean(size)) / self.size_rs_std(size)

    def zscore(self, set1, set2):
        # raw-score
        rs = self.hpo.get_hpo_set_eric_similarity(set1, set2)
        size = len(set1) * len(set2)
        if size == 0:
            return -999
        return (rs - self.size_rs_mean(size)) / self.size_rs_std(size)

    def predict_generator(self, hpo_list):
        hpo_list = list(filter(lambda x: x in self.hpo.hpo_set, hpo_list))
        for gene, gene_hpo_list in self.hpo.gene_hpo_series.iteritems():
            yield gene, self.zscore(hpo_list, gene_hpo_list)

    def predict(self, hpo_list, sort=False):
        # remove duplicated hpo
        hpo_list = list(set(hpo_list))
        result = pd.DataFrame(list(self.predict_generator(hpo_list)), columns=['gene_id', 'zscore'])
        if sort:
            result.sort_values(by='zscore', ascending=False, inplace=True)
            result.reset_index(drop=True, inplace=True)
        return result

    @staticmethod
    def norm(_value):
        if _value == 0:
            print('debug')
        MAX = 100
        MIN = -4
        _value = MAX if _value >= MAX else _value
        _value = MIN if _value <= MIN else _value
        if _value >= 0:
            return _value / (2 * MAX) + 0.5
        else:
            return (_value - MIN) / (2 * (0 - MIN))

    def predict_pkl(self, hpo_list, sort=True, add_parent=False, norm_value=False):
        """
        predict gene rank of specified hpo list
        :param hpo_list: hpo list
        :param sort: True to return sorted list
        :param add_parent: add parent hpo if hpo_list have empty gene list
        :param norm: normalize zscore to 0~1 to a new column score
        :return:
        """
        result_df = self.hpo.gene_hpo_sim_df[self.hpo.gene_hpo_sim_df.index.isin(set(hpo_list))]
        if len(result_df) == 0 and add_parent:
            hpo_list = list(
                it.chain.from_iterable(self.hpo.get_parent(x, valid_only=True) for x in hpo_list)
            )
            result_df = self.hpo.gene_hpo_sim_df[self.hpo.gene_hpo_sim_df.index.isin(hpo_list)]
        if len(result_df) == 0:
            warnings.warn(f'invalid hpo: {hpo_list}, zero will be return for all gene')
            result_df = result_df.transpose()
            result_df['zscore'] = 0
            result_df['rank'] = 0
            return result_df
        hpo_list = set(hpo_list)
        rs = result_df.sum()
        rs.name = 'rs'
        rs_df = pd.merge(rs, self.hpo.gene_hpo_size * len(hpo_list), left_index=True, right_index=True)
        rs_df['zscore'] = rs_df.apply(lambda x: self._zscore(x['rs'], x['size']), axis=1)
        rs_df.index.name = 'gene_id'
        if sort:
            rs_df.sort_values(by='zscore', ascending=False, inplace=True)
            rs_df['rank'] = range(1, len(rs_df) + 1)
        if norm_value:
            rs_df['score'] = rs_df['zscore'].map(self.norm)
        return rs_df

    def _hpo_zscore(self, _hpo_rs_seri, _hpo_list):
        _hpo_rs_seri[_hpo_list] = (_hpo_rs_seri[_hpo_list] - self.size_rs_mean(
            _hpo_rs_seri['size'])) / self.size_rs_std(len(_hpo_list) * _hpo_rs_seri['size'])
        return _hpo_rs_seri[_hpo_list]

    def predict_pkl_verbose(self, hpo_list, sort=True, return_valid_hpo=False, add_parent=False):
        if add_parent:
            hpo_map = dict((x, self.hpo.get_valid_hpo(x)) for x in set(hpo_list))
            hpo_list = it.chain.from_iterable(hpo_map.values())
        result_df = self.hpo.gene_hpo_sim_df[self.hpo.gene_hpo_sim_df.index.isin(hpo_list)]
        # if len(result_df) == 0 and add_parent:
        #     hpo_list = list(
        #         it.chain.from_iterable(self.hpo.get_parent(x, valid_only=True) for x in hpo_list)
        #     )
        #     result_df = self.hpo.gene_hpo_sim_df[self.hpo.gene_hpo_sim_df.index.isin(hpo_list)]
        new_hpo_list = result_df.index
        if len(result_df) == 0:
            raise ValueError(f'invalid hpo: {hpo_list}')

        rs_df = pd.merge(result_df.transpose(), self.hpo.gene_hpo_size, left_index=True, right_index=True)

        rs_df[new_hpo_list] = rs_df[new_hpo_list].sub(rs_df['size'].map(self.size_rs_mean), axis='index')\
            .div(rs_df['size'].map(lambda x: self.size_rs_std(x * len(new_hpo_list))), axis='index')
        rs_df['zscore'] = rs_df[new_hpo_list].sum(axis=1)
        rs_df.index.name = 'gene_id'
        rs_df = rs_df.merge(self.hpo.gene_df, how='left', left_index=True, right_index=True)
        if sort:
            rs_df.sort_values(by='zscore', ascending=False, inplace=True)
            rs_df['rank'] = range(1, len(rs_df) + 1)
        if return_valid_hpo:
            return list(new_hpo_list), rs_df
        return rs_df


if __name__ == '__main__':
    pass
    # _hpo = HPO()
    # psea = PSEA(_hpo)
    # hpo_li = ['HP:0008619', 'HP:0004322']
    # print('verbose', psea.predict_pkl_verbose(hpo_li).head())
    # print('ori', psea.predict_pkl(hpo_li).head())
