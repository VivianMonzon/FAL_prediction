import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--proteome_seq', required=True,
                    help='fasta file with protein sequences')
parser.add_argument('--out_feature', required=True,
                    help='output csv file containing ids and feature')
args = parser.parse_args()


def get_length(fh):
    id_len_dict = {}
    fh_in = open(fh)
    for name, seq in SimpleFastaParser(fh_in):
        length = len(seq)
        id_len_dict[name] = length
    df_id_len = pd.DataFrame(id_len_dict.items(), columns=['ID', 'length'])
    fh_in.close()
    return df_id_len


def treks(fh_in):
    df = pd.read_csv(fh_in, sep='\t')
    df['ID'] = df.seqid
    df = df[df.replength > 10]
    df['treks_07'] = 1
    df = df[['ID', 'treks_07']]
    df = df.fillna(0)
    return df


hmmsearch_cols = ['ID', 'tlen', 'qname', 'acc', 'qlen', 'Evalue', 'score', 'bias', 'cEvalue', 'iEvalue', 'Dscore', 'Dbias', 'from', 'to']


def anchor_search(fh_seq, hmm_anchor):
    fh = open(fh_seq)
    LPxTGs = []
    for name, seq in SimpleFastaParser(fh):
        seq1 = seq[-50:]
        sortase = re.compile('LP.T[G|A|N|D]')
        sortase2 = re.compile('NP.TG')
        sortase3 = re.compile('LP.GA')
        sortase4 = re.compile('LA.TG')
        sortase5 = re.compile('NPQTN')
        sortase6 = re.compile('IP.TG')
        if [m for m in re.finditer(sortase, seq1)]:
            LPxTGs.append(name)
        if [m for m in re.finditer(sortase2, seq1)]:
            LPxTGs.append(name)
        if [m for m in re.finditer(sortase3, seq1)]:
            LPxTGs.append(name)
        if [m for m in re.finditer(sortase4, seq1)]:
            LPxTGs.append(name)
        if [m for m in re.finditer(sortase5, seq1)]:
            LPxTGs.append(name)
        if [m for m in re.finditer(sortase6, seq1)]:
            LPxTGs.append(name)
        else:
            continue
    df = pd.read_csv(hmm_anchor, sep='\t', names=hmmsearch_cols)
    anchor_id = df.ID.tolist()
    ids_w_any_anchor = list(set(LPxTGs + anchor_id))
    df_anchor = pd.DataFrame({'ID': ids_w_any_anchor})
    df_anchor['Any_anchor'] = 1
    return df_anchor
    

def any_stalk(fh_in):
    df = pd.read_csv(fh_in, sep='\t', names=hmmsearch_cols)
    count_s = df['ID'].value_counts().rename_axis('ID').reset_index(name='counts')
    df_count_stalk = count_s.rename({'counts': 'Stalks'}, axis=1)
    return df_count_stalk


def any_adh(fh_in):
    df = pd.read_csv(fh_in, sep='\t', names=hmmsearch_cols)
    adhs = list(set(df.ID.tolist()))
    df_any_adh = pd.DataFrame({'ID': adhs})
    df_any_adh['Any_adh'] = 1
    return df_any_adh


def inmembrane(inmembrane_fh):
    df = pd.read_csv(inmembrane_fh, names=['ID', 'Prediction', 'Length', 'results', 'description'])
    cellwall_df = df[df['Prediction'] == 'PSE-Cellwall']
    cellwall_df['cellwall'] = 1
    PSE_df = cellwall_df[['ID', 'cellwall']]
    return PSE_df


def iupred(fh_in):
    df = pd.read_csv(fh_in)
    return df


def charge_hydro(fh_in):
    df = pd.read_csv(fh_in)
    return df


def aa_comp(fh_in):
    df = pd.read_csv(fh_in)
    return df


length = get_length(args.proteome_seq)

treks = treks('analysis/query_seq_treks.tsv')

df_any_anchor = anchor_search(args.proteome_seq,
                              'analysis/query_seq_anchor_dom_important.tbl')

continuous_stalks = any_stalk('analysis/query_seq_stalk_dom_important.tbl')
any_adhs = any_adh('analysis/query_seq_anchor_dom_important.tbl')

PSE = inmembrane('analysis/query_seq_PSE.csv')

iupred = iupred('analysis/query_seq_iupred_feature.csv')

hydro_charge = charge_hydro('analysis/query_seq_charge_hydro.csv')

aa_composition = aa_comp('analysis/query_seq_amino_acid_comp.csv')

df_all = length.merge(treks, how='left').merge(
    df_any_anchor, how='left').merge(continuous_stalks, how='left').merge(
        any_adhs, how='left').merge(PSE, how='left').merge(
            iupred, how='left').merge(hydro_charge, how='left').merge(
                aa_composition, how='left').drop_duplicates()
df_all = df_all.fillna(0)

df_all.to_csv(args.out_feature, index=False)
