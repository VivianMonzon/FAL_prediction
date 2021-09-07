import pandas as pd
import os


def preparations(seq, analysisfolder, resultsfolder):
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
    fh_seq = open('query_seq.fa', 'w')
    fh_seq.write('>query_sequence\n{}\n'.format(seq))
    fh_seq.close()


def hmmadh(hmmfolder):
    os.system('hmmsearch --cut_ga --domtblout analysis/query_seq_adh_dom.tbl {}/adh_dom_hmms.hmm query_seq.fa > analysis/query_seq_adh_dom.out'.format(hmmfolder))

    
class collect_features():
    
    def hmmsearch():
        os.system('hmmsearch --cut_ga --domtblout analysis/query_seq_adh_dom.tbl data/adh_dom_hmms.hmm query_seq.fa > analysis/query_seq_adh_dom.out')
        print('hmmsearch adh done')
        # python3.7 lib/make_out_readable.py --domtblout analysis/query_seq_adh_dom.tbl --out_readable analysis/query_seq_adh_dom_readable.tbl --out_important analysis/query_seq_adh_dom_important.tbl
