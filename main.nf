#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gigabiome/resisto-GEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Strain-specific genome-scale metabolic models from bacterial assemblies:
    AMRBench_QC gene prediction -> CarveMe -> species-reference biomass and gap-filling
    Github : https://github.com/gigabiome/resisto-GEM
----------------------------------------------------------------------------------
*/

include { RESISTO_GEM } from './workflows/resisto_gem'

workflow {
    RESISTO_GEM()
}
