import pandas as pd
import unittest
from io import StringIO
from hydro_charge_calculat import (proportion_hydro,
                                   proportion_charge,
                                   read_seq_summ_results)


sequence_type1 = ('MTLPLRRVLAAGGALALVAVLASPGVADPLPVDEPFSGTTFADPAWVLPGGDANSATLS'
                  'GDALRLTTPSPGFQVANATLDDPFRSDVAFTMDFDYGAYGG')
sequence_type2 = 'VSTVVTPPTPTPSASDPSPSSSA'
sequence_type3 = 'MTLPLRRVLAAXXVAVLASPGVAD'
sequence_type4 = 'UUUUUUUUUUUUUUUUUU'
sequence_type5 = 'ATTGCTTAAGCTAGGGCCAAT'

written_fasta1 = StringIO(
    '>sp|Q2G015|CLFA_STAA8 Clumping factor A OS=Staphylococcus '
    'aureus (strain NCTC 8325 / PS 47) OX=93061 GN=clfA PE=1\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
    '>sp|Q2G016|CLFB_STAA8 Clumping factor A OS=S. aureus '
    '(strain NCTC 8325 / PS 47) OX=93061 GN=clfB PE=1 SV=1\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENQQTQSDSASNESKS')
written_fasta2 = StringIO(
    '>Q2G015.1\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
    '>Q2G016.1\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n')
written_fasta3 = StringIO(
    '>Q2G015 Clumping factor A S. aureus\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
    '>Q2G016 Clumping factor B S. aureus\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS')
written_fasta4 = StringIO(
    '>Q2G015 Clumping factor A S. aureus\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS')
written_fasta5 = StringIO(
    '>Q2G015 Clumping factor A S. aureus\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS')
written_fasta6 = StringIO(
    '>Q2G015\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
    '>\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS')
written_fasta7 = StringIO(
    '>Q2G015\n'
    'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
    '>Q2G016\n'
    '')


class TestHydroCharge(unittest.TestCase):

    def test_hydrophobic_proportion(self):
        self.assertAlmostEqual(proportion_hydro(sequence_type1), 0.47)
        self.assertAlmostEqual(proportion_hydro(sequence_type2), 0.22)
        self.assertAlmostEqual(proportion_hydro(sequence_type3), 0.58)
        self.assertAlmostEqual(proportion_hydro(sequence_type4), 0)
        self.assertAlmostEqual(proportion_hydro(sequence_type5), 0.29)

    def test_charge_proportion(self):
        self.assertAlmostEqual(proportion_charge(sequence_type1), 0.15)
        self.assertAlmostEqual(proportion_charge(sequence_type2), 0.04)
        self.assertAlmostEqual(proportion_charge(sequence_type3), 0.12)
        self.assertAlmostEqual(proportion_charge(sequence_type4), 0)
        self.assertAlmostEqual(proportion_charge(sequence_type5), 0)
        
    def test_read_in_create_df(self):
        res1 = pd.DataFrame([['Q2G015', 0.34, 0.27], ['Q2G016', 0.32, 0.27]],
                            columns=['ID', 'Hydro_portion', 'Charge_portion'],
                            index=[0, 1])
        pd.testing.assert_frame_equal(read_seq_summ_results(written_fasta1),
                                      res1)
        res2 = pd.DataFrame([['Q2G015', 0.34, 0.27], ['Q2G016', 0.34, 0.27]],
                            columns=['ID', 'Hydro_portion', 'Charge_portion'],
                            index=[0, 1])
        pd.testing.assert_frame_equal(read_seq_summ_results(written_fasta2),
                                      res2)
        res3 = pd.DataFrame([['Q2G015', 0.34, 0.27], ['Q2G016', 0.34, 0.27]],
                            columns=['ID', 'Hydro_portion', 'Charge_portion'],
                            index=[0, 1])
        pd.testing.assert_frame_equal(read_seq_summ_results(written_fasta3),
                                      res3)
        res4 = pd.DataFrame([['Q2G015', 0.34, 0.27]],
                            columns=['ID', 'Hydro_portion', 'Charge_portion'],
                            index=[0])
        pd.testing.assert_frame_equal(read_seq_summ_results(written_fasta4),
                                      res4)
        res5 = pd.DataFrame([['Q2G015', 0.34, 0.27]],
                            columns=['ID', 'Hydro_portion', 'Charge_portion'],
                            index=[0])
        pd.testing.assert_frame_equal(read_seq_summ_results(written_fasta5),
                                      res5)
        self.assertRaises(ValueError, read_seq_summ_results, written_fasta6)
        self.assertRaises(ValueError, read_seq_summ_results, written_fasta7)


if __name__ == '__main__':
    unittest.main()
