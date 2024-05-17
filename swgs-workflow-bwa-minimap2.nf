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


log.info """\

S H A L L O W   W G S   C N - C A L L I N G  -  N F    v 2.1 
================================
genome   : $params.refFasta
reads    : $params.reads
results  : $params.outdir
"""


// Get reference genome
process DOWNLOAD_MM10 {
    errorStrategy 'retry', maxRetries: 2
    tag "Download Reference Genome mm10"
    
    output:
    path "reference_genome/mm10.fa.gz"

    script:
    """
    mkdir -p reference_genome
    curl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz > reference_genome/mm10.fa.gz
    """
}
process DOWNLOAD_HG19 {
    errorStrategy 'retry', maxRetries: 2
    tag "Download Reference Genome hg19"
    
    output:
    path "reference_genome/hg19.fa.gz"

    script:
    """
    mkdir -p reference_genome
    curl https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz > reference_genome/hg19.fa.gz
    """
}
process DOWNLOAD_HG38 {
    errorStrategy 'retry', maxRetries: 2
    tag "Download Reference Genome hg38"
    
    output:
    path "reference_genome/hg38.fa.gz"

    script:
    """
    mkdir -p reference_genome
    curl https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz > reference_genome/hg38.fa.gz
    """
}

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
    fastqc -o fastqc_${sample_id}_logs -f fastq $reads
    """
}

// Create a QC report
process MULTIQC1 {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'
 
    output:
    path "reports/pre_alignment_multiqc_report.html"
 
    script:
    """
    mkdir -p reports
    multiqc .
    mv multiqc_report.html reports/pre_alignment_multiqc_report.html
    """
}

// Index reference genome for bwa-mem2
process INDEX_BWAREF {
    tag "Index Reference Genome $genome"
    publishDir params.outdir, mode:'copy'

    input:
    file genome

    output:
    path "${genome.simpleName}_bwa_refgenome_index/${genome.simpleName}.*"

    script:
    """
    mkdir -p ${genome.simpleName}_bwa_refgenome_index
    cp ${genome} ${genome.simpleName}_bwa_refgenome_index/${genome}
    bwa-mem2 index ${genome} -p ${genome.simpleName}_bwa_refgenome_index/${genome}
    """
}

// Trim reads using Trimmomatic
process TRIMMOMATIC {
    tag "TRIMMOMATIC on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("trimmed/${sample_id}.{f,r}.paired.fastq"), 
    path("trimmed/${sample_id}.{f,r}.unpaired.fastq")

    script:
    if (params.pairedend)
        """
        mkdir trimmed
        trimmomatic PE -threads ${params.nthreads} -phred33 $reads \
            trimmed/${sample_id}.f.paired.fastq \
            trimmed/${sample_id}.f.unpaired.fastq \
            trimmed/${sample_id}.r.paired.fastq \
            trimmed/${sample_id}.r.unpaired.fastq \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 MAXINFO:100:0.5
        """
    else
        """
        mkdir trimmed
        trimmomatic SE -threads ${params.nthreads} -phred33 $reads \
            trimmed/${sample_id}.f.paired.fastq \
            trimmed/${sample_id}.f.unpaired.fastq \
            ILLUMINACLIP:TruSeq3-SE:2:30:10:2 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 MAXINFO:100:0.5
        """

}

// minimap2 alignment process
process MINIMAP2_ALIGNMENT {
    tag "MINIMAP2 ALIGNMENT on $sample_id"
    
    input:
    tuple val(sample_id), path(reads), path(unpaired_reads)
    path(reference)

    output:
    tuple val(sample_id), path("alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam")

    script:
    """
    mkdir -p "alignment_${sample_id}"
    minimap2 -a -x sr -Y -K 100M -t ${params.nthreads} ${reference} ${reads} | samtools view -hbS | \
            samtools sort -m 2G -@ ${params.nthreads} -o alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam
    samtools index -@ ${params.nthreads} alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam
    """
}

// minimap2 alignment process
process BWAMEM2_ALIGNMENT {
    tag "BWA-MEM 2 ALIGNMENT on $sample_id"
    
    input:
    tuple val(sample_id), path(reads), path(unpaired_reads)
    // path(reference)
    path(bwa_reference)

    output:
    tuple val(sample_id), path("alignment_${sample_id}/${sample_id}.se.bwa.sorted.bam")

    script:
    bwaList = bwa_reference.collect()
    fa_path = bwaList[4]                           // Access individual elements from the list
    unpairedList = unpaired_reads.collect()
    forward_unpaired = unpairedList[0]
    """
    mkdir -p "alignment_${sample_id}"
    bwa-mem2 mem -M -t ${params.nthreads} ${fa_path} ${reads} | samtools view -hbS | \
            samtools sort -m 2G -@ ${params.nthreads} -o alignment_${sample_id}/${sample_id}.bwa.sorted.bam
    samtools index -@ ${params.nthreads} alignment_${sample_id}/${sample_id}.bwa.sorted.bam

    if [ ${params.pairedend} ] && ![ ${params.useboth} ]
    then
            bwa-mem2 mem -M -t ${params.nthreads} ${fa_path} ${forward_unpaired} | samtools view -hbS | \
                samtools sort -m 2G -@ ${params.nthreads} -o alignment_${sample_id}/${sample_id}.up.bwa.sorted.bam
            samtools index -@ ${params.nthreads} alignment_${sample_id}/${sample_id}.up.bwa.sorted.bam

            # Grab forward strand from P.E. post_alignment
            samtools view -F 16 -o alignment_${sample_id}/${sample_id}.pe.f.bwa.sorted.bam \
                alignment_${sample_id}/${sample_id}.bwa.sorted.bam
            samtools merge -f -@ ${params.nthreads} alignment_${sample_id}/${sample_id}.se.merged.bam \
                alignment_${sample_id}/${sample_id}.up.bwa.sorted.bam \
                alignment_${sample_id}/${sample_id}.pe.f.bwa.sorted.bam
            # Re-Sort now merged BAM
            samtools sort -m 2G -@ ${params.nthreads} -o alignment_${sample_id}/${sample_id}.se.merged.sorted.bam \
                alignment_${sample_id}/${sample_id}.se.merged.bam
            samtools index -@ ${params.nthreads} alignment_${sample_id}/${sample_id}.se.merged.sorted.bam
            rm alignment_${sample_id}/${sample_id}.bwa.sorted.bam
            mv alignment_${sample_id}/${sample_id}.se.merged.sorted.bam alignment_${sample_id}/${sample_id}.bwa.sorted.bam
    fi
    """
}//.after(INDEX_BWAREF)

// Mark duplicate reads using picardtools
process MARKDUP {
    tag "MARKDUP on $sample_id"
    publishDir params.outdir, mode:'copy'
    
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("output_bams/${sample_id}.se.bwa.sorted.mkdup.bam"), 
    path("metrics_files/${sample_id}.marked_duplicates.metrics.txt")
    
    script:
    """
    mkdir -p "output_bams"
    mkdir -p "metrics_files"
    java "-Xmx16g" -jar /usr/picard/picard.jar MarkDuplicates \
      I=${bam} \
      O=output_bams/${sample_id}.se.bwa.sorted.mkdup.bam \
      M=metrics_files/${sample_id}.marked_duplicates.metrics.txt
    """
}

// Calculating a few samtools alignment metrics
process INDCOVFLAG {
    tag "Re-index and sort on $sample_id"
    
    input:
    tuple val(sample_id), path(bam), path(metrics)

    output:
    tuple val(sample_id), path("metrics_files/${sample_id}.se.coverageTable.tsv"), 
    path("metrics_files/${sample_id}.flagstat.txt")
    
    script:
    """
    mkdir -p metrics_files
    samtools index -@ ${params.nthreads} ${bam}
    samtools coverage ${bam} > metrics_files/${sample_id}.se.coverageTable.tsv
    samtools flagstat -@ ${params.nthreads} ${bam} > metrics_files/${sample_id}.flagstat.txt
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
    path "reports/final_multiqc_report.html"
 
    script:
    """
    mkdir -p reports
    multiqc .
    mv multiqc_report.html reports/final_multiqc_report.html
    """
}

// Call Copy-Numbers using QDNAseq
process CN_QDNA1 {
    publishDir params.outdir, mode:'copy'

    input:
    path(script)
    val(binsize)
    path(binannos_dir)
    val(bams)

    output:
    path("relative_cns/qdnaseq")

    script:
    """
    mkdir -p "relative_cns/qdnaseq"
    printf '%s\n' "${bams.join('\n')}" > bamfileslist.txt
    Rscript ${script} ${binsize} ${params.nthreads} relative_cns/qdnaseq bamfileslist.txt
    rm bamfileslist.txt
    """
}

workflow {

    // Make some channels for the needed params, reads files, and scripts
    binsizes_ch = Channel.from(params.binsizes)
    qdnaseq_script_ch = file("$projectDir/scripts/runQDNAseq.R")
    wisex_script_ch = file("$projectDir/scripts/runQDNAseq.R")
    
    if (params.pairedend) {
        reads_ch = Channel
            .fromFilePairs(params.reads, checkIfExists: true)
    } else {
        reads_ch = Channel
                .fromPath(params.reads, checkIfExists: true, type: 'file')
                .map { file -> tuple(file.simpleName, file) }
    } 

    // Get the reference genome
    if (params.genome == 'hg19') {
        reference_genome_ch = DOWNLOAD_HG19()
    } else if ( params.genome == 'hg38') {
        reference_genome_ch = DOWNLOAD_HG38()
    } else if (params.genome == 'mm10') {
        reference_genome_ch = DOWNLOAD_MM10()
    } else {
        fail "Reference genome is not recognized, please correct input string."
    }
    
    // reads_ch.view()
    fastqc1_ch = FASTQC1(reads_ch)
    // fastqc1_ch.view()
    MULTIQC1(fastqc1_ch.collect())

    if (params.pairedend) {
        trimmed_ch = TRIMMOMATIC(reads_ch)
    } else {
        trimmed_ch = TRIMMOMATIC(reads_ch)
    }

    if (params.aligner == 'bwamem2') {
        bwa_reference_genome_ch = INDEX_BWAREF(reference_genome_ch)
        bwa_reference_genome_ch.view()
        align_ch = BWAMEM2_ALIGNMENT(trimmed_ch, bwa_reference_genome_ch)
    } else if (params.aligner == 'minimap2') {
        align_ch = MINIMAP2_ALIGNMENT(trimmed_ch, reference_genome_ch)
    } else {
        fail "Aligner is not recognized, please correct input string. Options are 'bwamem2' or 'minimap2'."
    }

    markdup_ch = MARKDUP(align_ch)
    indcovflag_ch = INDCOVFLAG(markdup_ch)
    fastqc2_ch = FASTQC2(markdup_ch)

    if (params.runqdnaseq) {
        bamlist = markdup_ch.map { tuple -> tuple[1] }
        CN_QDNA1(qdnaseq_script_ch, binsizes_ch, params.binannos, bamlist.collect())
    }

    if (params.runwisex) {
        wisex_cns_ch = CN_WX1(wisex_script_ch, binsizes_ch)
    }

    // Extract the file path from each tuple in the channel
    metrics = markdup_ch.map { tuple -> tuple[2] }
    flagstats = indcovflag_ch.map { tuple -> tuple[2] }
    combined_ch = metrics.mix(flagstats, fastqc2_ch).collect()

    MULTIQC2(combined_ch)

}