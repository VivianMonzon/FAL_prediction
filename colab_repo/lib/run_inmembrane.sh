#!bin/bash                                                                                                                                                                                                            
FASTA=$1

main(){
    activate_py2
    run_inmem
    deactivate_py2
}

activate_py2(){
    source activate python2_env
}

run_inmem(){
    inmembrane_scan $FASTA
}

deactivate_py2(){
    conda deactivate
}

main
