#!/bin/bash

main(){
    predict_FAL
}

predict_FAL(){
    python3.7 ML_predict.py predict --fasta_seqs test_sequence.fasta \
	      --treks_dir /Users/vmonzon/Downloads \
	      --lipop_dir /Users/vmonzon/Downloads/software/LipoP1.0a \
	      --signalp_dir /Users/vmonzon/Downloads/software/signalp-5.0/bin \
	      --tmhmm_dir /Users/vmonzon/Downloads/software/tmhmm-2.0c/bin \
	      --analysisfolder analysistest \
	      --resultsfolder resultstest \
	      --jobname test
}

main
