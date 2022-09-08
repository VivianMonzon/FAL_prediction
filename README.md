# Discovery of bacterial Fibrillar Adhesins-like (FA-like) proteins
This repository comprises a random forest classification approach for the prediction of FA-like proteins in Gram positive bacteria.
Fibrillar adhesins are bacterial surface proteins, which play an important role in the bacterial pathogenesis [[1]](#1). The proteins of this novel defined protein class can bind
to proteins, carbohydrates or even ice crystalls and can enable the colonization of host cells. <br/>
## Requirements
The identification features used in this classification approach include the search for identical tandem repeats by a tool called T-REKS [[2]](#2) and the prediction of disordered regions by IUPred2A [[3]](#3). T-REKS is available to download [here](https://bioinfo.crbm.cnrs.fr/index.php?route=tools&tool=3). It requires muscle v3, which can be downloaded and installed by following the instructions [here](https://drive5.com/muscle/downloads_v3.htm) (Alternatively clustal can be used, whereby the lib/analyse_testing_data.sh has to be adapted accordingly). IUPRED is available to download [here](https://iupred2a.elte.hu/download_new).<br>

The other dependencies required for the code are listed in the 'requirements.txt' file. If a newer python version is used, please adapt the [hmmsearch file](lib/run_hmmsearch.sh) and [analysis run file](lib/analyse_testing_data.sh) accordingly. 

## Running the FA-like protein prediction approach
```
python3.7 ML_predict.py predict \
	      --fasta_seqs test_sequence.fasta \				# Give fasta file with sequence oof interest
              --treks_dir /Users/vmonzon/Downloads \				# Specify absolute path to directory containing the T-REKS tool
	      --iupred_dir /Users/vmonzon/Downloads/software \		    	# Specify path to directory containing the iupred script
              --analysisfolder ANALYSIS \                           	   	# Optional: if not given, will use 'analysis'
              --resultsfolder RESULTS \						# Optional: if not given, will use 'results'
              --jobname NAME                              	    		# Optional: if not given, will use fasta file name as jobname
```

## References
<a id="1">[1]</a>
Monzon, V., Lafita, A. & Bateman, A. 
Discovery of fibrillar adhesins across bacterial species.
BMC Genomics 22, 550 (2021). https://doi.org/10.1186/s12864-021-07586-2 <br>
<a id="2">[2]</a>
Julien Jorda, Andrey V. Kajava.
T-REKS: identification of Tandem REpeats in sequences with a K-meanS based algorithm.
Bioinformatics, Volume 25, Issue 20, 15 October 2009, Pages 2632–2638, https://doi.org/10.1093/bioinformatics/btp482 <br>
<a id="3">[3]</a>
Gábor Erdős, Zsuzsanna Dosztányi. 
Analyzing Protein Disorder with IUPred2A. 
Current Protocols in Bioinformatics 2020;70(1):e99. https://doi.org/10.1002/cpbi.99 <br>
