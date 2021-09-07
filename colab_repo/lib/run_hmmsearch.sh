# adhesive domains:
hmmsearch --cut_ga --domtblout analysis/query_seq_adh_dom.tbl ../resources/adh_dom_hmms.hmm $1 > analysis/query_seq_adh_dom.out
python3.7 lib/make_out_readable.py --domtblout analysis/query_seq_adh_dom.tbl --out_readable analysis/query_seq_adh_dom_readable.tbl --out_important analysis/query_seq_adh_dom_important.tbl
# stalk domains:
hmmsearch --cut_ga --domtblout analysis/query_seq_stalk_dom.tbl ../resources/stalk_dom_hmms.hmm $1 > analysis/query_seq_stalk_dom.out
python3.7 lib/make_out_readable.py --domtblout analysis/query_seq_stalk_dom.tbl --out_readable analysis/query_seq_stalk_dom_readable.tbl --out_important analysis/query_seq_stalk_dom_important.tbl
# anchor domains:
hmmsearch --cut_ga --domtblout analysis/query_seq_anchor_dom.tbl ../resources/anchor_dom_hmms.hmm $1 > analysis/query_seq_anchor_dom.out
python3.7 lib/make_out_readable.py --domtblout analysis/query_seq_anchor_dom.tbl --out_readable analysis/query_seq_anchor_dom_readable.tbl --out_important analysis/query_seq_anchor_dom_important.tbl
