#!/usr/bin/env python
from Bio.SeqIO.FastaIO import SimpleFastaParser
import statistics
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fasta_in', required=True)
parser.add_argument('--csv_out', required=True)
args = parser.parse_args()


hydrophobic_aa = {'A': 1, 'I': 1, 'L': 1, 'M': 1,
                  'F': 1, 'W': 1, 'Y': 1, 'V': 1}


def proportion_hydro(sequence):
    """Calculates proportion of hydrophobic aa of protein seq
    """
    hydro_portion = [hydrophobic_aa.get(aa, 0) for aa in sequence]
    hydro_portion = round(statistics.mean(hydro_portion), 2)
    return hydro_portion


def proportion_charge(sequence):
    """Calculates proportion of charged aa of protein seq
    """
    charged_aa = ['E', 'D', 'K', 'R']
    seq_charged_aa = [aa for aa in sequence if aa in charged_aa]
    charge_proportion = round(len(seq_charged_aa)/len(sequence), 2)
    return charge_proportion


def read_seq_summ_results(fh_input):
    """Reads in fasta file and splits name from seq to calculate hydro
    and charge proportion from seq. At the end it creates a dataframe.
    """
    for name, seq in SimpleFastaParser(fh_input):
        hydro = proportion_hydro(seq)
        charge = proportion_charge(seq)
    df = pd.DataFrame({'ID': [name], 'Hydro_portion': [hydro],
                       'Charge_portion': [charge]})
    return df


fh_in = open(args.fasta_in)
df_out = read_seq_summ_results(fh_in)
df_out.to_csv(args.csv_out, index=False)
