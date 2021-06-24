import pandas as pd
import unittest
from io import StringIO
from amino_acid_comp import (relative_entropy, amino_acid_comp)


class TestHydroCharge(unittest.TestCase):
    def setUp(self):
        """Testing different headers (written_fasta 1-3),
        files with just one protein (written_fasta4), 
        files with sequence split on severl lines (written_fasta5),
        error in sequence (written_fasta6/7), 
        and relative entropy of 0 (written_fasta8)
        """
        self.written_fasta1 = StringIO(
            '>sp|Q2G015|CLFA_STAA8 Clumping factor A OS=Staphylococcus '
            'aureus (strain NCTC 8325 / PS 47) OX=93061 GN=clfA PE=1\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
            '>sp|Q2G016|CLFB_STAA8 Clumping factor A OS=S. aureus '
            '(strain NCTC 8325 / PS 47) OX=93061 GN=clfB PE=1 SV=1\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENQQTQSDSASNESKS')
        self.written_fasta2 = StringIO(
            '>Q2G015.1\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
            '>Q2G016.1\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n')
        self.written_fasta3 = StringIO(
            '>Q2G015 Clumping factor A S. aureus\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
            '>Q2G016 Clumping factor B S. aureus\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS')
        self.written_fasta4 = StringIO(
            '>Q2G015 Clumping factor A S. aureus\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS')
        self.written_fasta5 = StringIO(
            '>Q2G015 Clumping factor A S. aureus\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS')
        self.written_fasta6 = StringIO(
            '>Q2G015\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
            '>\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS')
        self.written_fasta7 = StringIO(
            '>Q2G015\n'
            'MNMKKKEKHAIRKKSIGVASVLVGTLIGFGLLSSKEADASENSVTQSDSASNESKS\n'
            '>Q2G016\n'
            '')
        self.written_fasta8 = StringIO(
            '>Q2G015\n'
            'RHKDESTNQCGPAVILMFYWRHKDESTNQCGPAVILMFYWRHKDESTNQCGPAVILMFYW\n')

    def test_relative_entropy(self):
        res1 = pd.DataFrame([['Q2G015', 0.68], ['Q2G016', 0.62]],
                            columns=['ID', 'rel_entropy'],
                            index=[0, 1])
        pd.testing.assert_frame_equal(relative_entropy(self.written_fasta1),
                                      res1)
        res2_3 = pd.DataFrame([['Q2G015', 0.68], ['Q2G016', 0.68]],
                              columns=['ID', 'rel_entropy'],
                              index=[0, 1])
        pd.testing.assert_frame_equal(relative_entropy(self.written_fasta2),
                                      res2_3)
        pd.testing.assert_frame_equal(relative_entropy(self.written_fasta3),
                                      res2_3)
        res4_5 = pd.DataFrame([['Q2G015', 0.68]],
                              columns=['ID', 'rel_entropy'],
                              index=[0])
        pd.testing.assert_frame_equal(relative_entropy(self.written_fasta4),
                                      res4_5)
        pd.testing.assert_frame_equal(relative_entropy(self.written_fasta5),
                                      res4_5)
        self.assertRaises(ValueError, relative_entropy, self.written_fasta6)
        self.assertRaises(ValueError, relative_entropy, self.written_fasta7)
        res8 = pd.DataFrame([['Q2G015', 0.0]],
                            columns=['ID', 'rel_entropy'],
                            index=[0])
        pd.testing.assert_frame_equal(relative_entropy(self.written_fasta8),
                                      res8)
        
    def test_amino_acid_comp(self):
        res1 = pd.DataFrame([['Q2G015', 1.79, 1.79, 14.29, 3.57, 7.14,
                              19.64, 3.57, 5.36, 1.79, 0.0, 7.14,
                              0.0, 8.93, 7.14, 5.36, 7.14, 3.57,
                              1.79, 0.0, 0.0],
                             ['Q2G016', 1.79, 1.79, 14.29, 3.57, 7.14,
                              17.86, 3.57, 5.36, 5.36, 0.0, 7.14,
                              0.0, 8.93, 5.36, 5.36, 7.14, 3.57,
                              1.79, 0.0, 0.0]],
                            columns=['ID', 'R', 'H', 'K', 'D', 'E', 'S', 'T',
                                     'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I',
                                     'L', 'M', 'F', 'Y', 'W'],
                            index=[0, 1])
        pd.testing.assert_frame_equal(amino_acid_comp(self.written_fasta1),
                                      res1)
        res2 = pd.DataFrame([['Q2G015', 1.79, 1.79, 14.29, 3.57, 7.14,
                              19.64, 3.57, 5.36, 1.79, 0.0, 7.14,
                              0.0, 8.93, 7.14, 5.36, 7.14, 3.57,
                              1.79, 0.0, 0.0],
                             ['Q2G016', 1.79, 1.79, 14.29, 3.57, 7.14,
                              19.64, 3.57, 5.36, 1.79, 0.0, 7.14,
                              0.0, 8.93, 7.14, 5.36, 7.14, 3.57,
                              1.79, 0.0, 0.0]],
                            columns=['ID', 'R', 'H', 'K', 'D', 'E', 'S', 'T',
                                     'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I',
                                     'L', 'M', 'F', 'Y', 'W'],
                            index=[0, 1])
        pd.testing.assert_frame_equal(amino_acid_comp(self.written_fasta2),
                                      res2)
        pd.testing.assert_frame_equal(amino_acid_comp(self.written_fasta3),
                                      res2)
        res4_5 = pd.DataFrame([['Q2G015', 1.79, 1.79, 14.29, 3.57, 7.14,
                                19.64, 3.57, 5.36, 1.79, 0.0, 7.14,
                                0.0, 8.93, 7.14, 5.36, 7.14, 3.57,
                                1.79, 0.0, 0.0]],
                              columns=['ID', 'R', 'H', 'K', 'D', 'E', 'S', 'T',
                                       'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I',
                                       'L', 'M', 'F', 'Y', 'W'],
                              index=[0])
        pd.testing.assert_frame_equal(amino_acid_comp(self.written_fasta4),
                                      res4_5)
        pd.testing.assert_frame_equal(amino_acid_comp(self.written_fasta5),
                                      res4_5)
        self.assertRaises(ValueError, amino_acid_comp, self.written_fasta6)
        self.assertRaises(ValueError, amino_acid_comp, self.written_fasta7)
        res8 = pd.DataFrame([['Q2G015', 5.0, 5.0, 5.0, 5.0, 5.0,
                              5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                              5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                              5.0, 5.0, 5.0]],
                            columns=['ID', 'R', 'H', 'K', 'D', 'E', 'S', 'T',
                                     'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I',
                                     'L', 'M', 'F', 'Y', 'W'],
                            index=[0])
        pd.testing.assert_frame_equal(amino_acid_comp(self.written_fasta8),
                                      res8)


if __name__ == '__main__':
    unittest.main()
