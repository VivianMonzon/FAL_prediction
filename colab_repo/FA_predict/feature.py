import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
import re
from math import log2


def create_folder(analysisfolder, resultsfolder):
    if os.path.exists('{}'.format(analysisfolder)):
        print('{} folder exists'.format(analysisfolder))
    else:
        os.system('mkdir {}'.format(analysisfolder))
        print('{} folder created'.format(analysisfolder))
    if os.path.exists('{}'.format(resultsfolder)):
        print('{} folder exists'.format(resultsfolder))
    else:
        os.system('mkdir {}'.format(resultsfolder))
        print('{} folder created'.format(resultsfolder))
        

def query_seq_to_file(seq):
    fh_seq = open('query_seq.fa', 'w')
    fh_seq.write('>query_sequence\n{}\n'.format(seq))
    fh_seq.close()


hmmsearch_cols = ['ID', 'tlen', 'qname', 'acc', 'qlen', 'Evalue', 'score',
                  'bias', 'cEvalue', 'iEvalue', 'Dscore', 'Dbias', 'from',
                  'to']

    
class collect_features():

    # # maybe creating class object for each seq for aa_comp and etc. 
    # def __init__(self, name, seq):
    #     self.name = name
    #     self.seq = seq
        
    def adapt_hmmsearch(domtbl, output):
        fh_domtbl = open(domtbl, 'r')
        fh_out = open(output, 'w')
        for line in fh_domtbl:
            if line.startswith('#'):
                continue
            else:
                wo_gaps = re.sub(' +', ' ', line)
                w_tabs_list = wo_gaps.rsplit(' ')
                selection = [0, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 19, 20]
                out_list = [w_tabs_list[i] for i in selection]
                out_string = '\t'.join([str(x) for x in out_list])
                fh_out.write('{}\n'.format(out_string))
        fh_out.close()
        fh_domtbl.close()

    def hydro_charge(name, seq):
        hydro_aa = ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']
        charged_aa = ['E', 'D', 'K', 'R']
        hydro_fraction = [aa for aa in seq if aa in hydro_aa]
        hydro_fraction = round(len(hydro_fraction)/len(seq), 2)
        charge_fraction = [aa for aa in seq if aa in charged_aa]
        charge_fraction = round(len(charge_fraction)/len(seq), 2)
        df_hydro_charge = pd.DataFrame({'ID': [name],
                                        'Hydro_portion': [hydro_fraction],
                                        'Charge_portion': [charge_fraction]})
        return df_hydro_charge

    def aa_comp(name, seq):
        amino_acids = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G',
                       'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
        divergence = []
        aa_dict = {}
        for aa in amino_acids:
            P = seq.count(aa)/len(seq)
            if P != 0:
                D = P * log2(P/0.05)
                divergence.append(D)
            aa_dict[aa] = round(((seq.count(aa)/len(seq))*100), 2)
        entropy = round(sum(divergence), 2)
        entropy_df = pd.DataFrame({'ID': [name], 'rel_entropy': [entropy]})
        df_aa = pd.DataFrame([aa_dict])
        df_aa['ID'] = name
        df_aa_comp = entropy_df.merge(df_aa, on='ID', how='left')
        return df_aa_comp

    # def aa_comp_s(self):
    #     amino_acids = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G',
    #                    'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
    #     divergence = []
    #     aa_dict = {}
    #     for aa in amino_acids:
    #         P = self.seq.count(aa)/len(self.seq)
    #         if P != 0:
    #             D = P * log2(P/0.05)
    #             divergence.append(D)
    #         aa_dict[aa] = round(((self.seq.count(aa)/len(self.seq))*100), 2)
    #     entropy = round(sum(divergence), 2)
    #     entropy_df = pd.DataFrame({'ID': [self.name], 'rel_entropy': [entropy]})
    #     df_aa = pd.DataFrame([aa_dict])
    #     df_aa['ID'] = self.name
    #     df_aa_comp = entropy_df.merge(df_aa, on='ID', how='left')
    #     return df_aa_comp
    
    def adapt_iupred(iupred_results):
        df_iupred = pd.read_csv(iupred_results, sep='\t',
                                names=['SEQID', 'POS', 'RES', 'IUPRED2'],
                                comment='#')
        df_iupred['ID'] = df_iupred.SEQID.apply(
            lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
                x.split('.')[0] if '.' in x else (
                    x.split(' ')[0] if ' ' in x else x)))
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

    def treks(treks_results):
        df = pd.read_csv(treks_results, sep='\t')
        df['ID'] = df.seqid.apply(
            lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
                x.split('.')[0] if '.' in x else (
                    x.split(' ')[0] if ' ' in x else x)))
        df = df[df.replength > 10]
        df['treks_07'] = 1
        df = df[['ID', 'treks_07']]
        df = df.fillna(0)
        return df

    hmmsearch_cols = ['ID', 'tlen', 'qname', 'acc', 'qlen', 'Evalue', 'score',
                      'bias', 'cEvalue', 'iEvalue', 'Dscore', 'Dbias', 'from',
                      'to']

    def anchor_search(fh_seq, hmm_anchor):
        fh = open(fh_seq)
        LPxTGs = []
        for name, seq in SimpleFastaParser(fh):
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
        df['ID'] = df['ID'].apply(
            lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
                x.split('.')[0] if '.' in x else (
                    x.split(' ')[0] if ' ' in x else x)))
        anchor_id = df.ID.tolist()
        ids_w_any_anchor = list(set(LPxTGs + anchor_id))
        df_anchor = pd.DataFrame({'ID': ids_w_any_anchor})
        df_anchor['Any_anchor'] = 1
        return df_anchor

    def number_stalk(hmm_stalk):
        df = pd.read_csv(hmm_stalk, sep='\t', names=hmmsearch_cols)
        df['ID'] = df['ID'].apply(
            lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
                x.split('.')[0] if '.' in x else (
                    x.split(' ')[0] if ' ' in x else x)))
        count_s = df['ID'].value_counts().rename_axis('ID').reset_index(
            name='counts')
        df_count_stalk = count_s.rename({'counts': 'Stalks'}, axis=1)
        return df_count_stalk

    def any_adh(hmm_adh):
        df = pd.read_csv(hmm_adh, sep='\t', names=hmmsearch_cols)
        df['ID'] = df['ID'].apply(
            lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
                x.split('.')[0] if '.' in x else (
                    x.split(' ')[0] if ' ' in x else x)))
        adhs = list(set(df.ID.tolist()))
        df_any_adh = pd.DataFrame({'ID': adhs})
        df_any_adh['Any_adh'] = 1
        return df_any_adh

    def inmembrane(inmembrane_fh):
        df = pd.read_csv(inmembrane_fh, names=['ID', 'Prediction', 'Length',
                                               'results', 'description'])
        cellwall_df = df[df['Prediction'] == 'PSE-Cellwall']
        cellwall_df['cellwall'] = 1
        PSE_df = cellwall_df[['ID', 'cellwall']]
        PSE_df['ID'] = PSE_df['ID'].apply(
            lambda x: x.split('|')[1].split('|')[0] if '|' in x else (
                x.split('.')[0] if '.' in x else (
                    x.split(' ')[0] if ' ' in x else x)))
        return PSE_df
