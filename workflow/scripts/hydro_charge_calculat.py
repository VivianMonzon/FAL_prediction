#!/usr/bin/env python

"""
Calculates a set of properties from a protein sequence:
    - hydrophobicity (according to a particular scale)
    - total charge (at pH 7.4)
Author:
  Joao Rodrigues
  j.p.g.l.m.rodrigues@gmail.com
Adapted by:
    Vivian Monzon
"""

from Bio.SeqIO.FastaIO import SimpleFastaParser
import statistics
import pandas as pd

Fauchere_Pliska = {'A':  0.31, 'R': -1.01, 'N': -0.60,
                   'D': -0.77, 'C':  1.54, 'Q': -0.22,
                   'E': -0.64, 'G':  0.00, 'H':  0.13,
                   'I':  1.80, 'L':  1.70, 'K': -0.99,
                   'M':  1.23, 'F':  1.79, 'P':  0.72,
                   'S': -0.04, 'T':  0.26, 'W':  2.25,
                   'Y':  0.96, 'V':  1.22}
aa_charge = {'E': -1, 'D': -1, 'K': 1, 'R': 1}

hydrophobic_aa = {'A': 1, 'I': 1, 'L': 1, 'M': 1, 'F':1, 'W': 1, 'Y': 1, 'V': 1}


def hydrophobic_mean(sequence):
    hydro_for_aa = [Fauchere_Pliska.get(aa) for aa in sequence]
    hydro_for_aa = [x for x in hydro_for_aa if x is not None]
    return round(statistics.mean(hydro_for_aa), 2)


def proportion_hydro(sequence):
    hydro_portion = [hydrophobic_aa.get(aa, 0) for aa in sequence]
    return round(statistics.mean(hydro_portion), 2)


def calculate_charge(sequence, charge_dict=aa_charge):
    """Calculates the charge of the peptide sequence at pH 7.4                                                                                                                                              
    """
    sc_charges = [charge_dict.get(aa, 0) for aa in sequence]
    return round(statistics.mean(sc_charges), 2)


fh_in = open(snakemake.input[0])
ids = []
hydros = []
hydro_portion = []
charges = []
for name, seq in SimpleFastaParser(fh_in):
    if '|' in name:
        name = name.split('|')[1].split('|')[0]
    if ' ' in name:
        name = name.split(' ')[0]
    ids.append(name)
    hydros.append(hydrophobic_mean(seq))
    hydro_portion.append(proportion_hydro(seq))
    charges.append(calculate_charge(seq))


df = pd.DataFrame({'ID': ids, 'Hydrophobicity': hydros,
                   'Hydro_portion': hydro_portion, 'Charge': charges})
df = df[['ID', 'Hydro_portion', 'Charge']]
df.to_csv(snakemake.output[0], index=False)
