#!/bin/bash

fasta_seqs=$1
fasta_file=query_sequence.fa
analysis_folder=$2
results_folder=$3

main(){
    seq_string_to_file
    hmmersearch
    treks
    inmembrane
    hydro_charge 
    iupred
    amino_acid_composition
    combine_features
    # tidy_up
}

seq_string_to_file(){
    echo '>query_sequence' > $fasta_file
    echo $fasta_seqs >> $fasta_file
}


hmmersearch(){
    sh lib/run_hmmsearch.sh $fasta_file
}

treks(){
    java -Xmx4G -jar ../workflow/scripts/T-ReksHPC.jar -f $fasta_file -t $analysis_folder/query_seq_treks.tsv -a $analysis_folder/query_seq_treks.aln -m muscle -s 0.7 -S 5 -L 50 > $analysis_folder/query_seq_treks.log
}

inmembrane(){
    sh lib/run_inmembrane.sh $fasta_file
    mv query_sequence.csv $analysis_folder/query_seq_PSE.csv
    rm -r query_sequence/
}

hydro_charge(){
    python3.7 lib/hydro_charge_calculat.py --fasta_in $fasta_file --csv_out $analysis_folder/query_seq_charge_hydro.csv
}

iupred(){
    python3.7 lib/iupred2a.py $fasta_file "long" > $analysis_folder/query_seq_iupred.tab
    python3.7 lib/iupred_feature.py
}

amino_acid_composition(){
    python3.7 lib/amino_acid_comp.py --proteome_seq $fasta_file
}

combine_features(){
    python3.7 lib/combine_feature.py --proteome_seq $fasta_file --out_feature $analysis_folder/query_seq_feature.csv
}

tidy_up(){
    rm $fasta_file
    rm inmembrane.config
}

main
