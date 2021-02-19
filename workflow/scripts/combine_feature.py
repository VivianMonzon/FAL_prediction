import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re


def get_length(fh):
    id_len_dict = {}
    fh_in = open(fh)
    for name, seq in SimpleFastaParser(fh_in):
        if '|' in name:
            name = name.split('|')[1].split('|')[0]
        if ' ' in name:
            name = name.split(' ')[0]
        length = len(seq)
        id_len_dict[name] = length
    df_id_len = pd.DataFrame(id_len_dict.items(), columns=['ID', 'length'])
    fh_in.close()
    return df_id_len


def treks(fh_in):
    df = pd.read_csv(fh_in, sep='\t')
    df['ID'] = df.seqid.apply(lambda x: x.split('|')[1].split('|')[0] if '|' in x else x)
    df = df[df.replength > 10]
    df_high_psim = df[df.psim > 0.9]
    df_high_psim['treks_09'] = 1
    df['treks_07'] = 1
    df = df.merge(df_high_psim, on='ID', how='left')
    df = df[['ID', 'treks_07', 'treks_09']].drop_duplicates()
    df = df.fillna(0)
    return df


hmmsearch_cols = ['ID', 'tlen', 'qname', 'acc', 'qlen', 'Evalue', 'score', 'bias', 'cEvalue', 'iEvalue', 'Dscore', 'Dbias', 'from', 'to']


def motif_search(fh_in):
    fh = open(fh_in)
    LPxTGs = []
    for name, seq in SimpleFastaParser(fh):
        if '|' in name:
            name = name.split('|')[1].split('|')[0]
        if ' ' in name:
            name = name.split(' ')[0]
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
    return LPxTGs
    

def anchor_domains(fh_in, LPxTG_list):
    df = pd.read_csv(fh_in, sep='\t', names=hmmsearch_cols)
    df['ID'] = df['ID'].apply(lambda x: x.split('|')[1].split('|')[0] if '|' in x else x)
    LPxTGs = df[df.qname == 'Gram_pos_anchor'].ID.tolist()
    LPxTGs = list(set(LPxTGs))
    LPxTGs = list(set(LPxTGs + LPxTG_list))
    df_LPxTG = pd.DataFrame({'ID': LPxTGs})
    df_LPxTG['LPxTG'] = 1
    df_others = df[df.qname != 'Gram_pos_anchor']
    df_others['anchor'] = 1
    df_others = df_others[['ID', 'anchor']]
    return df_LPxTG, df_others


def list_to_df(l, feature):
    df = pd.DataFrame({'ID': l})
    df[feature] = 1
    return df


def any_stalk(fh_in):
    df = pd.read_csv(fh_in, sep='\t', names=hmmsearch_cols)
    df['ID'] = df['ID'].apply(lambda x: x.split('|')[1].split('|')[0] if '|' in x else x)
    count_s = df['ID'].value_counts().rename_axis('ID').reset_index(name='counts')
    df_count_stalk = count_s.rename({'counts': 'Stalks'}, axis=1)
    return df_count_stalk


def any_adh(fh_in):
    df = pd.read_csv(fh_in, sep='\t', names=hmmsearch_cols)
    df['ID'] = df['ID'].apply(lambda x: x.split('|')[1].split('|')[0] if '|' in x else x)
    adhs = list(set(df.ID.tolist()))
    df_any_adh = pd.DataFrame({'ID': adhs})
    df_any_adh['Any_adh'] =1
    return df_any_adh


def inmembrane(cellwall_in, membrane_in):
    cellwall_df = pd.read_csv(cellwall_in, names=['ID'])
    cellwall_df['cellwall'] = 1
    membrane_df = pd.read_csv(membrane_in, names=['ID'])
    membrane_df['membrane'] = 1
    PSE_df = cellwall_df.merge(membrane_df, on='ID', how='outer')
    PSE_df = PSE_df.fillna(0)
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


length = get_length(snakemake.input[0])

treks = treks(snakemake.input[1])

LPxTG_motifs = motif_search(snakemake.input[0])
LPxTG, other_anchors = anchor_domains(snakemake.input[2], LPxTG_motifs)

continuous_stalks = any_stalk(snakemake.input[3])
any_adhs = any_adh(snakemake.input[4])

PSE = inmembrane(snakemake.input[5], snakemake.input[6])

iupred = iupred(snakemake.input[7])

hydro_charge = charge_hydro(snakemake.input[8])

aa_composition = aa_comp(snakemake.input[9])

df_all = length.merge(treks, how='left').merge(other_anchors, how='left').merge(LPxTG, how='left').merge(continuous_stalks, how='left').merge(any_adhs, how='left').merge(PSE, how='left').merge(iupred, how='left').merge(hydro_charge, how='left').merge(aa_composition, how='left').drop_duplicates()
df_all = df_all.fillna(0)

df_all.to_csv(snakemake.output[0], index=False)
