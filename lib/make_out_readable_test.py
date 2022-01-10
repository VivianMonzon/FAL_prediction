import unittest
from make_out_readable import adapt_output
from io import StringIO
import os

written_tbl = StringIO(
    '#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord\n'
    '# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target\n'
    '#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------\n'
    'V6M7M1.1             -           1569 Collagen_bind        PF05737.13   131   3.9e-80  253.7  44.5   1   5   6.7e-27   1.3e-26   80.5   4.0     3   127   275   394   273   398 0.95 V6M7M1_9BACL Uncharacterized protein\n'
    'V6M7M1.1             -           1569 Collagen_bind        PF05737.13   131   3.9e-80  253.7  44.5   2   5   6.1e-13   1.2e-12   35.3   1.7     4   109   408   517   405   540 0.81 V6M7M1_9BACL Uncharacterized protein\n'
    'V6M7M1.1             -           1569 Collagen_bind        PF05737.13   131   3.9e-80  253.7  44.5   3   5   2.9e-11   5.7e-11   29.9   4.1    14   117   558   649   552   661 0.79 V6M7M1_9BACL Uncharacterized protein\n'
    'V6M7M1.1             -           1569 Collagen_bind        PF05737.13   131   3.9e-80  253.7  44.5   4   5   8.2e-30   1.6e-29   89.9   2.4     3   124   683   811   681   817 0.92 V6M7M1_9BACL Uncharacterized protein')


class TestMakeOutReadable(unittest.TestCase):
    def test_adapt_output(self):
        adapted_out = 'V6M7M1.1\t1569\tCollagen_bind\tPF05737.13\t131\t3.9e-80\t253.7\t44.5\t6.7e-27\t1.3e-26\t80.5\t4.0\t273\t398\n' \
            'V6M7M1.1\t1569\tCollagen_bind\tPF05737.13\t131\t3.9e-80\t253.7\t44.5\t6.1e-13\t1.2e-12\t35.3\t1.7\t405\t540\n' \
            'V6M7M1.1\t1569\tCollagen_bind\tPF05737.13\t131\t3.9e-80\t253.7\t44.5\t2.9e-11\t5.7e-11\t29.9\t4.1\t552\t661\n' \
            'V6M7M1.1\t1569\tCollagen_bind\tPF05737.13\t131\t3.9e-80\t253.7\t44.5\t8.2e-30\t1.6e-29\t89.9\t2.4\t681\t817\n'
        fh_tmp = open('tmp.tbl', 'w')
        adapt_output(written_tbl, fh_tmp)
        content = open('tmp.tbl').read()
        os.remove('tmp.tbl')
        self.assertEquals(content, adapted_out)

        
if __name__ == '__main__':
    unittest.main()
