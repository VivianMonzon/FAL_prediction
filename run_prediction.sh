#!/bin/bash

main(){
    predict_FAL
}

predict_FAL(){
    python3.7 ML_predict.py predict --fasta_seqs test_sequence.fasta \
	      --treks_dir /Users/vmonzon/Downloads \
	      --iupred_dir /Users/vmonzon/Downloads/software \
	      --analysisfolder analysistest2 \
	      --resultsfolder resultstest2 \
	      --jobname test2
}

main
