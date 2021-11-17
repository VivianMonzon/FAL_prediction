import pandas as pd
import argparse


def iupred(fh_in):
    df_iupred = pd.read_csv(fh_in, sep='\t',
                            names=['SEQID', 'POS', 'RES', 'IUPRED2'],
                            comment='#')
    df_iupred['ID'] = df_iupred.SEQID.apply(
        lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
            x.split('.')[0] if '.' in x else x))
    protein_ids = list(set(df_iupred.ID.tolist()))
    fraction_disorder = {}
    for x in protein_ids:
        df_one = df_iupred[df_iupred.ID == x]
        df_one_len = df_one.shape[0]
        perc_dis = round(((df_one[df_one.IUPRED2 > 0.5].shape[0]) / (df_one_len)) * 100, 2)
        fraction_disorder[x] = perc_dis
    df_iupred = pd.DataFrame(fraction_disorder.items(),
                             columns=['ID', 'frac_disordered'])
    df_iupred = df_iupred.drop_duplicates()
    return df_iupred


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fh_in', required=True, help='iupred results')
    parser.add_argument('--fh_out', required=True, help='df output')
    args = parser.parse_args()
    iupred_out = open(args.fh_in)
    df_out = iupred(iupred_out)
    df_out.to_csv(args.fh_out, index=False)
