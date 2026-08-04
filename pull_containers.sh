#!/usr/bin/env bash
set -euo pipefail

IMAGE_DIR=/projects/marralab/rcorbett_prj/singularity_cache
APPTAINER_CACHE=/projects/marralab/scratch/rcorbett/apptainer_cache
APPTAINER_TMP=/projects/marralab/scratch/rcorbett/apptainer_tmp

mkdir -p "$IMAGE_DIR" "$APPTAINER_CACHE" "$APPTAINER_TMP"

export APPTAINER_CACHEDIR="$APPTAINER_CACHE"
export APPTAINER_TMPDIR="$APPTAINER_TMP"

# Compatibility with installations still invoked as `singularity`
export SINGULARITY_CACHEDIR="$APPTAINER_CACHE"
export SINGULARITY_TMPDIR="$APPTAINER_TMP"

pull_image() {
    local output_name="$1"
    local source_uri="$2"
    local output_path="${IMAGE_DIR}/${output_name}"

    echo
    echo "============================================================"
    echo "Image:  $source_uri"
    echo "Output: $output_path"
    echo "============================================================"

    # Remove only abandoned partial files for this image.
    rm -f "${output_path}.pulling."*

    if [[ -s "$output_path" ]]; then
        if singularity inspect "$output_path" >/dev/null 2>&1; then
            echo "Already present and valid; skipping."
            return 0
        fi

        echo "Existing image failed inspection; moving it aside."
        mv "$output_path" "${output_path}.invalid.$(date +%Y%m%d_%H%M%S)"
    fi

    singularity pull "$output_path" "$source_uri"
    singularity inspect "$output_path" >/dev/null

    echo "Pull and inspection succeeded."
}

pull_image \
    clinicalgenomics-bwa-mem2-2.2.1.img \
    docker://clinicalgenomics/bwa-mem2:2.2.1

pull_image \
    huntsmanlab-azure_download-latest.img \
    docker://huntsmanlab/azure_download:latest

pull_image \
    huntsmanlab-ichorcna-latest.img \
    docker://huntsmanlab/ichorcna:latest

pull_image \
    huntsmanlab-utanos-latest.img \
    docker://huntsmanlab/utanos:latest

pull_image \
    huntsmanlab-qdnaseq-latest.img \
    docker://huntsmanlab/qdnaseq:latest

pull_image \
    sofvdvel-wisecondorx-0.1.img \
    docker://sofvdvel/wisecondorx:0.1

pull_image \
    staphb-fastqc-latest.img \
    docker://staphb/fastqc:latest

pull_image \
    broadinstitute-picard-latest.img \
    docker://broadinstitute/picard:latest

pull_image \
    quay.io-biocontainers-multiqc-1.3--py35_2.img \
    docker://quay.io/biocontainers/multiqc:1.3--py35_2

pull_image \
    staphb-samtools-1.19.img \
    docker://staphb/samtools:1.19

pull_image \
    staphb-trimmomatic-latest.img \
    docker://staphb/trimmomatic:latest

echo
echo "All requested images were pulled and successfully inspected."

echo
echo "Image checksums:"
sha256sum "$IMAGE_DIR"/*.img


