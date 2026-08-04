#!/usr/bin/env nextflow

include { SWGS_QC_PIPELINE } from './workflows/swgs_qc'

workflow {
    log.info """\

    sWGS PROCESSING PIPELINE  v 3.0
    ================================
    output directory = ${params.output_directory}
    reference genome = ${params.reference_genome_version}
    """

    WGS_QC_PIPELINE()
}
