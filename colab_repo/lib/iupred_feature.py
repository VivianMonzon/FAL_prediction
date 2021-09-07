import pandas as pd


def iupred(fh_in):
    df_iupred = pd.read_csv(fh_in, sep='\t',
                            names=['SEQID', 'POS', 'RES', 'IUPRED2'],
                            comment='#')
    df_iupred['ID'] = df_iupred.SEQID
    protein_ids = list(set(df_iupred.ID.tolist()))
    fraction_disorder = {}
    for x in protein_ids:
        df_one = df_iupred[df_iupred.ID == x]
        df_one_len = df_one.shape[0]
        perc_dis = round(((
            df_one[df_one.IUPRED2 > 0.5].shape[0]) / (df_one_len)) * 100, 2)
        fraction_disorder[x] = perc_dis
    iupred_frac = pd.DataFrame(fraction_disorder.items(),
                               columns=['ID', 'frac_disordered'])
    return iupred_frac


df_out = iupred('analysis/query_seq_iupred.tab')
df_out.to_csv('analysis/query_seq_iupred_feature.csv', index=False)
