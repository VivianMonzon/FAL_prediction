# Discovery of bacterial Fibrillar Adhesins-like (FA-like) proteins
This repository contains a snakemake pipeline to predict FA-like proteins in Gram positive bacteria using a RandomForest classification approach.
Fibrillar adhesins are bacterial surface proteins, which play an important role in the bacterial pathogenesis [[1]](#1). The proteins of this novel defined protein class can bind
to proteins, carbohydrates or even ice crystalls and can enable the colonization of host cells. <br/>
To predict FA-like proteins, indicate your sequence file in the config/config.yaml file and run snakemake using the provided conda environments (optional: number of cores to use):<br/>
```bash
snakemake --use-conda --cores 4
```
Alternatively, the google colab notebook can be used (colab_repo/ML_FA_prediction.ipynb).

## References
<a id="1">[1]</a>
Monzon, V., Lafita, A. & Bateman, A. 
Discovery of fibrillar adhesins across bacterial species.
BMC Genomics 22, 550 (2021). https://doi.org/10.1186/s12864-021-07586-2