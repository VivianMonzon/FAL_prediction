#!bin/bash

FASTA=$1
lipop_dir=$2
signalp_dir=$3
tmhmm_dir=$4
analysis_dir=$5

main(){
    write_config
    activate_py2
    run_inmem
    deactivate_py2
}

write_config(){
    echo "{\n  'fasta': '',\n  'csv': '',\n  'out_dir': '',\n  'protocol': 'gram_pos',\n\n " \
    "'signalp4_bin': '$signalp_dir/signalp',\n  'lipop1_bin': '$lipop_dir/LipoP',\n  'tmhmm_bin': '$tmhmm_dir/tmhmm',\n  'memsat3_bin': 'runmemsat',\n  'helix_programs': ['tmhmm'],\n  'terminal_exposed_loop_min': 50,\n  'internal_exposed_loop_min': 100,\n\n " \
    "'hmmsearch3_bin': 'hmmsearch',\n  'hmm_evalue_max': 0.1,\n  'hmm_score_min': 10,\n\n " \
    "'barrel_programs': ['tmbetadisc-rbf'],\n  'bomp_clearly_cutoff': 3,\n  'bomp_maybe_cutoff': 1,\n  'tmbetadisc_rbf_method': 'aadp'}" > inmembrane.config
}

activate_py2(){
    ENVS=$(conda env list | grep 'inmembrane_env' | cut -d' ' -f1)
    if [[ $ENVS = *"inmembrane_env"* ]]; then
    	echo 'Environment exists'
	continue
    else
    	echo "inmembrane_env doesnot exists - will create it"
    	conda env create -f env/inmembrane_env.yml
    fi;
    source activate inmembrane_env
}

run_inmem(){
    inmembrane_scan $FASTA &> $analysis_dir/inmembrane.log
}

deactivate_py2(){
    conda deactivate
}

main