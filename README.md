# ML_predict_snakemake_pipeline
Snakemake pipeline to predict FA-like proteins based on RandomForest algorithm
Test it with the test_sequence.fasta file (if other file change it in config/config.yaml)
Run: $ snakemake -np (n: dry run) --use-conda (to use conda environment specified) --cores 8 (number of cores)
