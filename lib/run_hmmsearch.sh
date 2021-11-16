# adhesive domains:
hmmsearch --cut_ga --domtblout $3/${2}_adh_dom.tbl database/adh_dom_hmms.hmm $1 > $3/${2}_adh_dom.out
python3.7 lib/make_out_readable.py --fh_in $3/${2}_adh_dom.tbl --fh_out $3/${2}_adh_dom_adapted.tbl
# stalk domains:
hmmsearch --cut_ga --domtblout $3/${2}_stalk_dom.tbl database/stalk_dom_hmms.hmm $1 > $3/${2}_stalk_dom.out
python3.7 lib/make_out_readable.py --fh_in $3/${2}_stalk_dom.tbl --fh_out $3/${2}_stalk_dom_adapted.tbl
# anchor domains:
hmmsearch --cut_ga --domtblout $3/${2}_anchor_dom.tbl database/anchor_dom_hmms.hmm $1 > $3/${2}_anchor_dom.out
python3.7 lib/make_out_readable.py --fh_in $3/${2}_anchor_dom.tbl --fh_out $3/${2}_anchor_dom_adapted.tbl
