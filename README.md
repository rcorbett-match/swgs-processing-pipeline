# sWGS Processing Pipeline using Nextflow and Singularity
### Author: Maxwell Douglas

------------------------------------------------------------------------

## Introduction

### Background/Rationale

Shallow whole genome sequencing (sWGS) can be used to detect copy number (CN) aberrations, detect CN-Signatures (1), detect Homologous Recombination Defficiency (HRD) (2), and even create CN-Signatures (1). This is a popular sequencing type used for neo-natal diagnostics and studying cancer. One of the major benefits is the reduced cost when compared with Whole Genome Sequencing. With just 15 million reads or so usually targeted this assay costs a fraction of Whole Genome Sequencing. This nextflow workflow serves as a pre-processing pipeline for sWGS data and can be the first step before doing some interesting downstream analyses. Several recent papers have also specifically made use of a sample-prep protocol that uses single-ended sequencing (50bp) from FFPE tissue (1). This single ended protocol isn't quite as typical in this day and age of NGS where paired-end 150bp is usually the default. This workflow therefore serves a unique purpose in being specifically with that use-case in mind.  
In the guide that follows below we will assume the reader has some basic familiarity with using the command-line (ex. bash shell), installing software using the commandline, and using git for version control/interacting with github.

### Order of operations

This pipeline will cover the following pre-processing steps:

1. Download a reference genome (curl) - this step isn't really pre-processing, but a reference genome is needed in order to do alignment.
2. Sequencing quality assessment of the input reads. (FastQC)
3. Aggregation of these inital QC reports for each sample into a single report. (MultiQC)
4. Alignment of the reads (minimap2)
5. Duplicate read identification and marking (Picard tools)
6. Assess read coverage and alignment statistics (Samtools)
7. Quality assessment of the now aligned reads (FastQC)
8. Aggregation of this second round of QC reports for each sample into a single report. (MultiQC)

### Visual representation of pipeline execution (DAG)

![](flowchart.png)

### Directories/Files in this repository

-   `test_data` : Directory containing some test data that can be used with this pipeline. This test data is not published though yet and needs to be kept in a private repo for now.
-   `nextflow.config`: Config file for Nextflow, contain all Docker container info and setting needed for using Singularity for this pipeline
-   `workflow.nf`: The main pipeline file that contain the instructions for running the pipeline through Nextflow.
-   `misc_scripts`: Old, stub, or under-development scripts. This directory can be safely ignored.
-   `.gitignore`: What files/directories to ignore when developing and using git as the version control system.
-   `README.md`: This file!

## Usage

### Software Installation

In order to use this pipeline, a user must have `git` , `nextflow`, and `singularity` all properly installed with appropriate permissions for the user. This pipeline does NOT require sudo access.  

If any of the above is not yet installed please follow their respective guides:  
For [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)  
For [nextflow](https://www.nextflow.io/docs/latest/getstarted.html)  
For [singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)  
  
Next, since this pipeline is in a private repository, ensure that you have created the appropriate personal access tokens on github.  
[This guide](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens) provides appropriate instructions  

Last, install this pipeline by cloning the repo.  
```         
git clone https://github.com/genomaxx/swgs-processing-pipeline.git
```

### How to run this workflow

From the commandline (assuming you have just cloned this git repo), navigate into the newly created directory.  
`cd swgs-processing-pipeline`  

To run the pipeline on the sample data, simply now execute the following command:  
`nextflow run swgs-workflow-se.nf -resume`  

To use your own data, open the `nextflow.config` file and replace the path in the 5th line with the path to your own data.  
So, for example:  
`reads = "$baseDir/test_data/chr21_raw/**.{fastq,fq,fastq.gz,fq.gz}"`  
may become...  
`reads = "$baseDir/../yournewdatadirectory/somefileprefix**.{fastq,fq,fastq.gz,fq.gz}"`  
Then save your changes, and execute the pipeline using the `nextflow run` command as seen above.

Look for the `results` directory (the output) to appear in this same run directory after pipeline execution.  

## Input
This pipeline expects single-ended raw short-read sequencing as input. (ex. from Illumina) The reads are expected in `fastq` formated files with any one of the following extensions: `XXX.fastq` | `XXX.fastq.gz` | `XXX.fq` | `XXX.fq.gz`  
The pipeline expects one file per sample. The fastq format is a common data standard who's details can be found [on wikipedia](https://en.wikipedia.org/wiki/FASTQ_format).
A brief outline of that formatting is copied below for convenience.

    A FASTQ file has four line-separated fields per sequence:
    Field 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
    Field 2 is the raw sequence letters.
    Field 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
    Field 4 encodes the quality values for the sequence in Field 2, and must contain the same number of symbols as letters in the sequence.
    A FASTQ file containing a single sequence might look like this:
    
    @SEQ_ID
    GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
    +
    !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65

## Output
A folder named `results` will be created with the execution of this pipeline.  
It's contents will include:  

1. Two QC reports - one pre-alignment:  
   `pre_alignment_multiqc_report.html`  
   and one post-alignment:  
   `final_multiqc_report.html`  
   both `.html` files.  
   Here is an image of what those files should look like when opened up in a browser.  
![](multiqc_report.png)
   
   For detailed insutructions on how to interpret the plots generated in the MultiQC reports there are many great resources online.  
   [Here is one.](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html)
   
3. A folder named `output_bams` containing the aligned reads. One file per sample.  
   These files will be in the `bam` format and more details on their formatting can be found on the [wikipedia page](https://en.wikipedia.org/wiki/Binary_Alignment_Map).
   Also contained within this folder are files named like such: `XXX.marked_duplicates.metrics.txt`, they contain metrics about the number of duplicated reads in each sample.
   There will be of these marked_duplicates files per input sample.

## Extras

### Troubleshooting
In order for this pipeline to work as expected, the software installations must be done properly. 
This is especially important for singularity and nextflow (they have many dependencies and are compplex pieces of software).  
For example, on one of the testing machines used in developing this pipeline the following commands were needed for singularity to run as expected:  
`source /cvmfs/soft.computecanada.ca/config/profile/bash.sh`
`module load apptainer`
`module load nextflow`

### List of containers used by singularity in this workflow
`docker://curlimages/curl:latest`
`docker://staphb/fastqc:latest`  
`docker://quay.io/biocontainers/multiqc:1.3--py35_2`  
'docker://niemasd/minimap2_samtools:latest'    
`docker://broadinstitute/picard:latest`  
`docker://staphb/samtools:latest`  
`docker://staphb/fastqc:latest`  
`docker://quay.io/biocontainers/multiqc:1.3--py35_2`  

### References

1. Macintyre, G., Goranova, T. E., De Silva, D., Ennis, D., Piskorz, A. M., Eldridge, M., Sie, D., Lewsley, L. A., Hanif, A., Wilson, C., Dowson, S., Glasspool, R. M., Lockley, M., Brockbank, E., Montes, A., Walther, A., Sundar, S., Edmondson, R., Hall, G. D., Clamp, A., … Brenton, J. D. (2018). Copy number signatures and mutational processes in ovarian carcinoma. Nature genetics, 50(9), 1262–1270. https://doi.org/10.1038/s41588-018-0179-8
2. Callens, C., Rodrigues, M., Briaux, A. et al. Shallow whole genome sequencing approach to detect Homologous Recombination Deficiency in the PAOLA-1/ENGOT-OV25 phase-III trial. Oncogene 42, 3556–3563 (2023). https://doi.org/10.1038/s41388-023-02839-8
