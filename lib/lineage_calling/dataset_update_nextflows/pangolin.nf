// pangolin wrapper script, everything is through conda, does nothing

nextflow.enable.dsl=2

process update_caller_dataset {
    cpus 1

    shell:
    '''
    #! /usr/bin/env bash

    exit 0;
    '''
    }


workflow {
    update_caller_dataset()
}