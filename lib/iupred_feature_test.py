import pandas as pd
import unittest
from io import StringIO
from iupred_feature import iupred

written_fasta1 = StringIO(
    '# IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding\n'
    '# Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi\n'
    '# Nucleic Acids Research 2018;46(W1):W329-W337.\n'
    '#\n'
    '# Prediction type: long\n'
    '# Prediction output\n'
    '# SEQID\tPOS\tRES\tIUPRED2\n'
    'V6M7M1.1\t1\tM\t0.2609\n'
    'V6M7M1.1\t2\tW\t0.2364\n'
    'V6M7M1.1\t3\tE\t0.2164\n')
written_fasta2 = StringIO(
    '# IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding\n'
    '# Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi\n'
    '# Nucleic Acids Research 2018;46(W1):W329-W337.\n'
    '#\n'
    '# Prediction type: long\n'
    '# Prediction output\n'
    '# SEQID\tPOS\tRES\tIUPRED2\n'
    'V6M7M1.1\t1\tM\t0.2609\n'
    'V6M7M1.1\t2\tW\t0.2364\n'
    'V6M7M1.1\t3\tE\t0.2164\n'
    'V6M7M2.1\t1\tM\t0.2609\n'
    'V6M7M2.1\t2\tW\t0.2364\n'
    'V6M7M2.1\t3\tE\t0.7164\n'
    'V6M7M2.1\t4\tG\t0.7609\n'
    'V6M7M2.1\t5\tA\t0.8364\n'
    'V6M7M2.1\t6\tG\t0.9164\n')


class TestIupredAdapt(unittest.TestCase):

    def test_iupred(self):
        res1 = pd.DataFrame([['V6M7M1', 0.0]],
                            columns=['ID', 'frac_disordered'],
                            index=[0])
        pd.testing.assert_frame_equal(iupred(written_fasta1),
                                      res1)
        res2 = pd.DataFrame([['V6M7M2', 66.67],
                             ['V6M7M1', 0.0]],
                            columns=['ID', 'frac_disordered'],
                            index=[0, 1])
        pd.testing.assert_frame_equal(iupred(written_fasta2),
                                      res2)


if __name__ == '__main__':
    unittest.main()
