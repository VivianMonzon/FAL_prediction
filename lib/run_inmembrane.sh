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
    echo "{
'fasta': '',
'csv': '',
'out_dir': '',
'protocol': 'gram_pos',
                                                                                                                                                                                                                                                                              
'signalp4_bin': '$signalp_dir/signalp',
'lipop1_bin': '$lipop_dir/LipoP',
'tmhmm_bin': '$tmhmm_dir/tmhmm',
'memsat3_bin': 'runmemsat',
'helix_programs': ['tmhmm'],
'terminal_exposed_loop_min': 50,
'internal_exposed_loop_min': 100,

'hmmsearch3_bin': 'hmmsearch',
'hmm_evalue_max': 0.1,
'hmm_score_min': 10,

'barrel_programs': ['tmbetadisc-rbf'],
'bomp_clearly_cutoff': 3,
'bomp_maybe_cutoff': 1,
'tmbetadisc_rbf_method': 'aadp'
}" > inmembrane.config    
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
