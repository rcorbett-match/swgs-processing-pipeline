#!/usr/bin/env nextflow

/*
 * Copyright (c) 2020, Maxwell Douglas
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 * 
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

/* 
 * 'Shallow WGS - NF' - A Nextflow pipeline for processing shallow WGS data
 * 
 * This pipeline reproduces the analysis pipeline described in " Targeted 
 * and shallow whole genome sequencing identifies therapeutic opportunities 
 * in p53abn endometrial cancers"
 * 
 * Maxwell Douglas 
 */

params.reads = "$baseDir/test_data/chr1_raw/**.{fastq,fq,fastq.gz,fq.gz}"
params.outdir = "$baseDir/results"
params.refFasta = "$baseDir/../reference_genome/hg19/GRCh37-lite.fa"
params.bwaThreads = 8
params.haplotypeCallerMappingQual = 20
params.haplotypeCallerMinBaseQual = 20

// Pre-alignment QC
process FASTQC {
    tag "FASTQC on $sample_id"
    publishDir params.outdir

    input:
    tuple val(sample_id), path(reads)
 
    output:
    path "fastqc_${sample_id}_logs"
 
    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc "$reads" -o fastqc_${sample_id}_logs
    """
}



// Create a QC report
process MULTIQC {
    publishDir params.outdir, mode:'copy'

    container "docker://quay.io/biocontainers/multiqc:1.3--py35_2"

    input:
    path '*'
 
    output:
    path "multiqc_report.html"
 
    script:
    """
    multiqc .
    """
}


workflow {
    log.info """\

S H A L L O W   W G S   C N - C A L L I N G  -  N F    v 2.1
================================
genome   : $params.refFasta
reads    : $params.reads
results  : $params.outdir
"""

    reads_ch = Channel
                .fromPath( params.reads, checkIfExists: true, type: 'file' )
                .map { file -> tuple(file.simpleName, file) }

    fastqc_ch=FASTQC(reads_ch)
    // align_ch=ALIGNMENT(reads_ch)
    MULTIQC(fastqc_ch.collect())
}
