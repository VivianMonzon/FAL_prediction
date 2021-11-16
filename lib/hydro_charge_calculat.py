#!/usr/bin/env python
from Bio.SeqIO.FastaIO import SimpleFastaParser
import statistics
import pandas as pd
import argparse

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
    ids = []
    hydro_portion = []
    charge_proportion = []
    for name, seq in SimpleFastaParser(fh_input):
        if '|' in name:
            name = name.split('|')[1].split('|')[0]
        if '.' in name:
            name = name.split('.')[0]
        if ' ' in name:
            name = name.split(' ')[0]
        if name in ids:
            raise ValueError('File contains proteins with same ID!')
        if name == '':
            raise ValueError('Header couldnot select ID!')
        ids.append(name)
        if not seq:
            raise ValueError('Sequence is empty')
        hydro_portion.append(proportion_hydro(seq))
        charge_proportion.append(proportion_charge(seq))
    df = pd.DataFrame({'ID': ids, 'Hydro_portion': hydro_portion,
                       'Charge_portion': charge_proportion})
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fh_in', required=True, help='Fasta file')
    parser.add_argument('--fh_out', required=True, help='df output')
    args = parser.parse_args()
    seq_fh = open(args.fh_in)
    df_out = read_seq_summ_results(seq_fh)
    df_out.to_csv(args.fh_out, index=False)

    
# fh_in = open(snakemake.input[0])
# df_out.to_csv(snakemake.output[0], index=False)

