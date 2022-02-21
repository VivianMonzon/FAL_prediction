#!/bin/bash

fasta_seqs=$1
analysis_folder=$2
results_folder=$3
jobname=$4
treks=$5
iupred=$9

main(){
    hmmersearch
    treks
    hydro_charge
    iupred                                                                                                                                                                                                 
    amino_acid_composition
    combine_features
}

hmmersearch(){
    echo 'Domain hmmsearch'
    sh lib/run_hmmsearch.sh $fasta_seqs $jobname $analysis_folder
}

treks(){
    echo 'T-REKS tandem repeat search'
    java -Xmx4G -jar $treks/T-ReksHPC.jar -f $fasta_seqs -t $analysis_folder/${jobname}_treks.tsv \
	 -a $analysis_folder/${jobname}_treks.aln -m muscle -s 0.7 -S 5 -L 50 > $analysis_folder/${jobname}_treks.log
}

hydro_charge(){
    echo 'Calculation hydrophobic and charged aa'
    python3.7 lib/hydro_charge_calculat.py --fh_in $fasta_seqs \
	      --fh_out $analysis_folder/${jobname}_charge_hydro.csv
}

iupred(){
    echo 'Iupred predictions and feature collection'
    python3.7 $iupred/iupred2a.py $fasta_seqs "long" > $analysis_folder/${jobname}_iupred.tab
    python3.7 lib/iupred_feature.py --fh_in $analysis_folder/${jobname}_iupred.tab \
	      --fh_out $analysis_folder/${jobname}_iupred_feature.csv
}

amino_acid_composition(){
    echo 'AA composition calculation'
    python3.7 lib/amino_acid_comp.py --fh_in $fasta_seqs --fh_out $analysis_folder/${jobname}_amino_acid_comp.csv
}

combine_features(){
    echo 'Combining features'
    python3.7 lib/combine_feature.py --fasta_fh $fasta_seqs \
	      --treks_fh $analysis_folder/${jobname}_treks.tsv \
	      --anchor_fh $analysis_folder/${jobname}_anchor_dom_adapted.tbl \
	      --stalk_fh $analysis_folder/${jobname}_stalk_dom_adapted.tbl \
	      --adh_fh $analysis_folder/${jobname}_adh_dom_adapted.tbl \
	      --iupred_fh $analysis_folder/${jobname}_iupred_feature.csv \
	      --hydro_charge_fh $analysis_folder/${jobname}_charge_hydro.csv \
	      --aa_comp_fh $analysis_folder/${jobname}_amino_acid_comp.csv \
	      --out_fh $analysis_folder/${jobname}_feature.csv
}

main
