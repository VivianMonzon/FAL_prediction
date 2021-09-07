import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
from math import log2

parser = argparse.ArgumentParser()
parser.add_argument('--proteome_seq', required=True, help='fasta file with protein sequences')
args = parser.parse_args()


amino_acids = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']


def relative_entropy(fasta_file):
    fh = open(fasta_file)
    dict_id_diverg = {}
    for name, seq in SimpleFastaParser(fh):
        divergence = []
        for aa in amino_acids:
            P = seq.count(aa)/len(seq)
            if P != 0:
                D = P * log2(P/0.05)
                divergence.append(D)
        entropy = sum(divergence)
        dict_id_diverg[name] = round(entropy, 2)
    df_entropy = pd.DataFrame(dict_id_diverg.items(),
                              columns=['ID', 'rel_entropy'])
    return df_entropy


def amino_acid_comp(fh_fasta):
    fh = open(fh_fasta)
    columns = ['ID'] + amino_acids
    df_all = pd.DataFrame(columns=columns, dtype=object)
    for name, seq in SimpleFastaParser(fh):
        aa_dict = {}
        for aa in amino_acids:
            aa_dict[aa] = round(((seq.count(aa)/len(seq))*100), 2)
        aa_df = pd.DataFrame([aa_dict])
        aa_df['ID'] = name
        df_all = df_all.append(aa_df, sort=False)
    return df_all


rel_entropy = relative_entropy(args.proteome_seq)
aa_composition = amino_acid_comp(args.proteome_seq)
merged_df = rel_entropy.merge(aa_composition, on='ID', how='outer')
merged_df.to_csv('analysis/query_seq_amino_acid_comp.csv', index=False)
