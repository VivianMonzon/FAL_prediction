import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from math import log2

all_amino_acids = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']


def relative_entropy(fasta_file):
    fh = open(fasta_file)
    dict_id_diverg = {}
    for name, seq in SimpleFastaParser(fh):
        if '|' in name:
            name = name.split('|')[1].split('|')[0]
        if ' ' in name:
            name = name.split(' ')[0]
        divergence = []
        for aa in all_amino_acids:
            P = seq.count(aa)/len(seq)
            if P != 0:
                D = P * log2(P/0.05)
                divergence.append(D)
        entropy = sum(divergence)
        dict_id_diverg[name] = round(entropy, 2)
    df_entropy = pd.DataFrame(dict_id_diverg.items(), columns=['ID', 'rel_entropy'])
    return df_entropy


selected_amino_acids = ['T', 'L']


def amino_acid_comp(fh_fasta):
    fh = open(fh_fasta)
    df_all = pd.DataFrame()
    for name, seq in SimpleFastaParser(fh):
        if '|' in name:
            name = name.split('|')[1].split('|')[0]
        if ' ' in name:
            name = name.split(' ')[0]
        aa_dict = {}
        for aa in selected_amino_acids:
            aa_dict[aa] = round(((seq.count(aa)/len(seq))*100), 2)
        aa_df = pd.DataFrame([aa_dict])
        aa_df['ID'] = name
        df_all = df_all.append(aa_df, sort=False)
    return df_all


rel_entropy = relative_entropy(snakemake.input[0])
aa_composition = amino_acid_comp(snakemake.input[0])
merged_df = rel_entropy.merge(aa_composition, on='ID', how='outer')
merged_df.to_csv(snakemake.output[0], index=False)
