import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re
import argparse


def get_length(fh_in):
    id_len_dict = {}
    for name, seq in SimpleFastaParser(fh_in):
        if '|' in name:
            name = name.split('|')[1].split('|')[0]
        if '.' in name:
            name = name.split('.')[0]
        if ' ' in name:
            name = name.split(' ')[0]
        length = len(seq)
        id_len_dict[name] = length
    df_id_len = pd.DataFrame(id_len_dict.items(), columns=['ID', 'length'])
    return df_id_len


def treks(fh_in):
    df = pd.read_csv(fh_in, sep='\t')
    df['ID'] = df.seqid.apply(lambda x: x.split('|')[1].split('|')[0]
                              if '|' in x else (
                                      x.split('.')[0] if '.' in x else (
                                          x.split(' ')[0] if ' ' in x else x)))
    df = df[df.replength > 10]
    df['treks_07'] = 1
    df = df[['ID', 'treks_07']].drop_duplicates()
    df = df.fillna(0)
    return df


hmmsearch_cols = ['ID', 'tlen', 'qname', 'acc', 'qlen', 'Evalue', 'score', 'bias', 'cEvalue', 'iEvalue', 'Dscore', 'Dbias', 'from', 'to']


def anchor_search(fh_seq, hmm_anchor):
    LPxTGs = []
    for name, seq in SimpleFastaParser(fh_seq):
        if '|' in name:
            name = name.split('|')[1].split('|')[0]
        if ' ' in name:
            name = name.split(' ')[0]
        if '.' in name:
            name = name.split('.')[0]
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
    df['ID'] = df['ID'].apply(lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
        x.split('.')[0] if '.' in x else (x.split(' ')[0] if ' ' in x else x)))
    anchor_id = df.ID.tolist()
    ids_w_any_anchor = list(set(LPxTGs + anchor_id))
    df_anchor = pd.DataFrame({'ID': ids_w_any_anchor})
    df_anchor['Any_anchor'] = 1
    return df_anchor


def list_to_df(l, feature):
    df = pd.DataFrame({'ID': l})
    df[feature] = 1
    return df


def any_stalk(fh_in):
    df = pd.read_csv(fh_in, sep='\t', names=hmmsearch_cols)
    df['ID'] = df['ID'].apply(lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
        x.split('.')[0] if '.' in x else (x.split(' ')[0] if ' ' in x else x)))
    count_s = df['ID'].value_counts().rename_axis('ID').reset_index(
        name='counts')
    df_count_stalk = count_s.rename({'counts': 'Stalks'}, axis=1)
    return df_count_stalk


def any_adh(fh_in):
    df = pd.read_csv(fh_in, sep='\t', names=hmmsearch_cols)
    df['ID'] = df['ID'].apply(lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
        x.split('.')[0] if '.' in x else (x.split(' ')[0] if ' ' in x else x)))
    adhs = list(set(df.ID.tolist()))
    df_any_adh = pd.DataFrame({'ID': adhs})
    df_any_adh['Any_adh'] = 1
    return df_any_adh


def inmembrane(inmembrane_fh):
    df = pd.read_csv(inmembrane_fh, names=['Name', 'Prediction', 'Length',
                                           'results', 'description'])
    cellwall_df = df[df['Prediction'] == 'PSE-Cellwall']
    cellwall_df['cellwall'] = 1
    cellwall_df['ID'] = cellwall_df['Name'].apply(
        lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
            x.split('.')[0] if '.' in x else (
                x.split(' ')[0] if ' ' in x else x)))
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


# if __name__ == '__main__':
#     fh_seq = open(snakemake.input[0])
#     length = get_length(fh_seq)
#     fh_seq.close()
#     treks = treks(snakemake.input[1])
#     fh_seq = open(snakemake.input[0])
#     any_anchor_df = anchor_search(fh_seq, snakemake.input[2])
#     continuous_stalks = any_stalk(snakemake.input[3])
#     fh_seq.close()
#     any_adhs = any_adh(snakemake.input[4])
#     PSE = inmembrane(snakemake.input[5])
#     iupred = iupred(snakemake.input[6])
#     hydro_charge = charge_hydro(snakemake.input[7])
#     aa_composition = aa_comp(snakemake.input[8])
#     df_all = length.merge(treks, how='left').merge(any_anchor_df, how='left').merge(continuous_stalks, how='left').merge(any_adhs, how='left').merge(PSE, how='left').merge(iupred, how='left').merge(hydro_charge, how='left').merge(aa_composition, how='left').drop_duplicates()
#     df_all = df_all.fillna(0)
#     df_all.to_csv(snakemake.output[0], index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_fh', required=True,
                        help='Fasta file')
    parser.add_argument('--treks_fh', required=True)
    parser.add_argument('--anchor_fh', required=True)
    parser.add_argument('--stalk_fh', required=True)
    parser.add_argument('--adh_fh', required=True)
    parser.add_argument('--pse_fh', required=True)
    parser.add_argument('--iupred_fh', required=True)
    parser.add_argument('--hydro_charge_fh', required=True)
    parser.add_argument('--aa_comp_fh', required=True)
    parser.add_argument('--out_fh', required=True)
    args = parser.parse_args()
    fh_seq = open(args.fasta_fh)
    length = get_length(fh_seq)
    fh_seq.close()
    treks = treks(args.treks_fh)
    fh_seq = open(args.fasta_fh)
    any_anchor_df = anchor_search(fh_seq, args.anchor_fh)
    continuous_stalks = any_stalk(args.stalk_fh)
    fh_seq.close()
    any_adhs = any_adh(args.adh_fh)
    PSE = inmembrane(args.pse_fh)
    iupred = iupred(args.iupred_fh)
    hydro_charge = charge_hydro(args.hydro_charge_fh)
    aa_composition = aa_comp(args.aa_comp_fh)
    df_all = length.merge(treks, how='left').merge(
        any_anchor_df, how='left').merge(continuous_stalks, how='left').merge(
            any_adhs, how='left').merge(PSE, how='left').merge(
                iupred, how='left').merge(hydro_charge, how='left').merge(
                    aa_composition, how='left').drop_duplicates()
    df_all = df_all.fillna(0)
    df_all.to_csv(args.out_fh, index=False)
