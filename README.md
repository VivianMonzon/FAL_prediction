# Discovery of bacterial Fibrillar Adhesins-like (FA-like) proteins
This repository contains a snakemake pipeline to predict FA-like proteins in Gram positive bacteria using a RandomForest classification approach.
Fibrillar adhesins are bacterial surface proteins, which play an important role in the bacterial pathogenesis [[1]](#1). The proteins of this novel defined protein class can bind
to proteins, carbohydrates or even ice crystalls and can enable the colonization of host cells. <br/>
## Requirements
The identification features used in this classification approach include the search for identical tandem repeats by a tool called T-REKS [[2]](#2) and the search for cell surface anchors by the pipeline inmembrane [[3]](#3). T-REKS is available to download [here](https://bioinfo.crbm.cnrs.fr/index.php?route=tools&tool=3). Conda is required for our classification approach to use the environment provided to run inmembrane. Additionally, Inmembrane is build on the following tools:<br>

<ul>
<li>TMHMM 2.0</li>
<ul>
<li> Download: https://services.healthtech.dtu.dk/cgi-bin/sw_request </li>
<li> Instructions: https://ssbio.readthedocs.io/en/latest/instructions/tmhmm.html </li>
</ul>
<li>SignalP 4.1</li>
<ul>
<li> Download: https://services.healthtech.dtu.dk/service.php?SignalP-4.1 </li>
</ul>
<li>LipoP 1.0</li>
<ul>
<li> Download: https://services.healthtech.dtu.dk/service.php?LipoP-1.0 </li>
</ul>
</ul>


## Running the FA-like protein prediction approach
```
python3.7 ML_predict.py predict \
	      --fasta_seqs test_sequence.fasta \				# Give fasta file with sequence oof interest
              --treks_dir /Users/vmonzon/Downloads \				# Specify absolute path to T-REKS directory
              --lipop_dir /Users/vmonzon/Downloads/software/LipoP1.0a \         # Specify absolute path to LipoP directory
              --signalp_dir /Users/vmonzon/Downloads/software/signalp-5.0/bin \ # Specify absolute path to SingalP directory
              --tmhmm_dir /Users/vmonzon/Downloads/software/tmhmm-2.0c/bin \    # Specify absolute path to TMHMM directory
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
Perry, A.J., Ho, B.K.
Inmembrane, a bioinformatic workflow for annotation of bacterial cell-surface proteomes.
Source Code Biol Med 8, 9 (2013). https://doi.org/10.1186/1751-0473-8-9