#!/usr/bin/env python 
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from math import log2
import argparse


amino_acids = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']


def relative_entropy(fasta_file):
    """investigating the sequence composition bias by calculating the rel. Entropy.
    A relative entropy of 0 would mean that all amino acids are eqally present.
    """
    dict_id_diverg = {}
    for name, seq in SimpleFastaParser(fasta_file):
        if ' ' in name:
            name = name.split(' ')[0]
        if '|' in name:
            name = name.split('|')[1].split('|')[0]
        if '.' in name:
            name = name.split('.')[0]
        if name == '':
            raise ValueError('Header couldnot select ID!')
        if not seq:
            raise ValueError('Sequence is empty')
        divergence = []
        for aa in amino_acids:
            P = seq.count(aa)/len(seq)
            if P != 0:
                D = P * log2(P/0.05)
                divergence.append(D)
        entropy = sum(divergence)
        dict_id_diverg[name] = round(entropy, 2)
    df_entropy = pd.DataFrame(dict_id_diverg.items(), columns=['ID',
                                                               'rel_entropy'])
    return df_entropy


def amino_acid_comp(fh_fasta):
    """Calculating for each amino acid the proportion it is present in the 
    regarding sequence.
    """
    columns = ['ID'] + amino_acids
    df_all = pd.DataFrame(columns=columns, dtype=object)
    for name, seq in SimpleFastaParser(fh_fasta):
        if ' ' in name:
            name = name.split(' ')[0]
        if '|' in name:
            name = name.split('|')[1].split('|')[0]
        if '.' in name:
            name = name.split('.')[0]
        if name == '':
            raise ValueError('Header couldnot select ID!')
        if not seq:
            raise ValueError('Sequence is empty')
        aa_dict = {}
        for aa in amino_acids:
            aa_dict[aa] = round(((seq.count(aa)/len(seq))*100), 2)
        aa_df = pd.DataFrame([aa_dict])
        aa_df['ID'] = name
        df_all = df_all.append(aa_df, sort=False)
    df_all = df_all.reset_index(drop=True)
    return df_all


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fh_in', required=True, help='Fasta file')
    parser.add_argument('--fh_out', required=True, help='df output')
    args = parser.parse_args()
    fh_seq = open(args.fh_in)
    df_rel_E = relative_entropy(fh_seq)
    fh_seq.close()
    fh_seq = open(args.fh_in)
    df_aa = amino_acid_comp(fh_seq)
    merged_df = df_rel_E.merge(df_aa, on='ID', how='outer')
    merged_df.to_csv(args.fh_out, index=False)
    

# rel_entropy = relative_entropy(snakemake.input[0])
# aa_composition = amino_acid_comp(snakemake.input[0])
# merged_df = rel_entropy.merge(aa_composition, on='ID', how='outer')
# merged_df.to_csv(snakemake.output[0], index=False)
