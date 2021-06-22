# ML_predict_snakemake_pipeline
Snakemake pipeline to predict Fibrillar adhesins in bacterial proteins using a RandomForest classification approach.<br/>
Fibrillar adhesins are bacterial surface proteins, which play an important role in the bacterial pathogenesis. The proteins of this novel defined protein class can bind
to proteins, carbohydrates or even ice crystalls and can enable the colonization of host cells. They are characterised by their protein architecture,
particularly their repeating protein domains, also called stalk domains, with which they fold into a stalk to project the adhesive domain closer to the binding target.<br/>
To test the pipeline a protein sequence example file is provided ('test_sequence.fasta'). To use your own sequence, change the name of the to analyse fasta file in the config/config.yaml file.<br/>
For the analysis run snakemake using the given conda enviroments as indicated in the Snakefile, e.g.:<br/>
```bash
snakemake --use-conda --cores 8
```