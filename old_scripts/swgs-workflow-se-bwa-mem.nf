#!/usr/bin/env nextflow

/*
 * Copyright (c) 2023, Maxwell Douglas
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


// Pre-alignment QC
process FASTQC1 {
    tag "FASTQC on $sample_id"

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
process MULTIQC1 {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'
 
    output:
    path "pre_alignment_multiqc_report.html"
 
    script:
    """
    multiqc .
    mv multiqc_report.html pre_alignment_multiqc_report.html
    """
}

// Alignment step
process ALIGNMENT {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam") 
    // path("alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam.bai")

    script:
    """
    mkdir -p "alignment_${sample_id}"
    bwa-mem2 mem -M -t ${params.bwaThreads} ${params.refFasta} \
        ${reads} > alignment_${sample_id}/${sample_id}.se.bwa.sam
    samtools view -h -b -S alignment_${sample_id}/${sample_id}.se.bwa.sam > alignment_${sample_id}/${sample_id}.se.bwa.bam
    samtools sort -m 2G -@ ${params.bwaThreads} -o alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam alignment_${sample_id}/${sample_id}.se.bwa.bam
    samtools index -@ ${params.bwaThreads} alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam
    rm alignment_${sample_id}/${sample_id}.se.bwa.sam
    rm alignment_${sample_id}/${sample_id}.se.bwa.bam
    """
}

// MARKDUP process
process MARKDUP {
    tag "MARKDUP on $sample_id"
    publishDir params.outdir, mode:'copy'
    
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("markdup/${sample_id}.se.bwa.sorted.mkdup.bam"), 
    path("markdup/${sample_id}.marked_duplicates.metrics.txt")
    
    script:
    """
    mkdir -p "markdup"
    java "-Xmx16g" -jar /usr/picard/picard.jar MarkDuplicates \
      I=${bam} \
      O=markdup/${sample_id}.se.bwa.sorted.mkdup.bam \
      M=markdup/${sample_id}.marked_duplicates.metrics.txt
    """
}

// MARKDUP process
process INDCOVFLAG {
    tag "Re-index and sort on $sample_id"
    
    input:
    tuple val(sample_id), path(bam), path(metrics)

    output:
    tuple val(sample_id), path("indcovflag/${sample_id}.se.coverageTable.tsv"), 
    path("indcovflag/${sample_id}.flagstat.txt")
    
    script:
    """
    mkdir -p indcovflag
    samtools index -@ 24 ${bam}
    samtools coverage ${bam} > indcovflag/${sample_id}.se.coverageTable.tsv
    samtools flagstat -@ ${params.bwaThreads} ${bam} > indcovflag/${sample_id}.flagstat.txt
    """
}

// Post-alignment QC
process FASTQC2 {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(bam), path(metrics)
 
    output:
    path "post_alignment_fastqc_logs_${sample_id}"
 
    script:
    """
    mkdir -p post_alignment_fastqc_logs_${sample_id}
    fastqc ${bam} -o post_alignment_fastqc_logs_${sample_id}
    """
}

// Create a QC report
process MULTIQC2 {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'
 
    output:
    path "final_multiqc_report.html"
 
    script:
    """
    multiqc .
    mv multiqc_report.html final_multiqc_report.html
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

    fastqc1_ch = FASTQC1(reads_ch)
    MULTIQC1(fastqc1_ch.collect())

    align_ch = ALIGNMENT(reads_ch)
    markdup_ch = MARKDUP(align_ch)
    indcovflag_ch = INDCOVFLAG(markdup_ch)
    fastqc2_ch = FASTQC2(markdup_ch)

    // Extract the file path from each tuple in the channel
    metrics = markdup_ch.map { tuple -> tuple[2] }
    flagstats = indcovflag_ch.map { tuple -> tuple[2] }
    combined_ch = metrics.mix(flagstats, fastqc2_ch).collect()

    MULTIQC2(combined_ch)

    combined_ch.subscribe { data ->
        // Process or print the combined data
        println "Combined data: $data"
    }

    test = fastqc1_ch.collect()
    test.subscribe { data ->
        // Process or print the combined data
        println "Combined data: $data"
    }


    // MULTIQC2(fastqc_ch.collect())

}

