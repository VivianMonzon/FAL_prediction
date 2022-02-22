import pandas as pd
import unittest
from io import StringIO
from combine_feature import (get_length, treks, anchor_search)


class TestCombine(unittest.TestCase):

    def test_length(self):
        written_fasta = StringIO(
            '>Q2G015.1\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
            '>Q2G016.1\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGF\n')
        res1 = pd.DataFrame([['Q2G015', 56], ['Q2G016', 29]],
                            columns=['ID', 'length'],
                            index=[0, 1])
        pd.testing.assert_frame_equal(get_length(written_fasta),
                                      res1)
        # self.assertRaises(ValueError, relative_entropy, self.written_fasta6)
        
    def test_treks(self):
        written_tsv = StringIO('seqid\trepnumber\treplength\tstart\tend\t'
                               'psim\ttotlength\n'
                               'Q2G015\t3\t30\t200\t350\t0.83\t250\n')
        written_tsv2 = StringIO('seqid\trepnumber\treplength\tstart\tend\t'
                                'psim\ttotlength\n'
                                'Q2G016\t15\t5\t300\t350\t0.89\t125\n')
        written_tsv3 = StringIO('seqid\trepnumber\treplength\tstart\tend\t'
                                'psim\ttotlength\n')
        res1 = pd.DataFrame([['Q2G015', 1]], columns=['ID', 'treks_07'],
                            index=[0])
        pd.testing.assert_frame_equal(treks(written_tsv), res1)
        res_empty = pd.DataFrame({}, columns=['ID', 'treks_07'],
                                 index=[], dtype='int64')
        assert treks(written_tsv2).shape == res_empty.shape
        assert treks(written_tsv3).shape == res_empty.shape

    def test_motif_search(self):
        written_fasta = StringIO(
            '>Q2G015.1\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
            '>Q2G016.1\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFLPTTG\n'
            '>Q2G017.1\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSLPSNESKS\n')
        written_tsv = StringIO(
            'Q2G015.1\tGram_pos_anchor\tPF00746.21\t43\t5e-09\t22.3\t0.0\t'
            '1e-08\t21.3\t0.0\t1566\n')
        res = pd.DataFrame([['Q2G015', 1], ['Q2G016', 1]],
                           columns=['ID', 'Any_anchor'], index=[0, 1])
        assert anchor_search(written_fasta, written_tsv).shape == res.shape


if __name__ == '__main__':
    unittest.main()
