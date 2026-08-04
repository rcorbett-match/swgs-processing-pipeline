# Nextflow 26 Behavioral Equivalence Checklist

This checklist covers the compatibility migration to Nextflow `26.04.6`.

| Area | Status | Notes |
| --- | --- | --- |
| Processes | Unchanged | No process names were changed. |
| Process commands | Unchanged except one compile-blocking name fix | Commands were not refactored. `QDNA_BINS` changed an undefined interpolation from `qd_bwavgbed` to documented `params.qd_bwgavgbed`. |
| Process inputs and outputs | Unchanged | No process `input:` or `output:` blocks were edited. |
| Workflow invocation order | Preserved | `main.nf` now invokes the included `SWGS_QC_PIPELINE`; internal call order was preserved. Old script top-level statements were moved into their existing workflow blocks before existing calls. |
| Channel construction and operators | Preserved | Existing channel operators were not changed. WisecondorX process calls now expand existing emitted outputs by index to satisfy strict arity checking. |
| Parameters and defaults | Unchanged | `config/params.yaml` and config defaults were not edited. |
| Containers and environments | Unchanged | FastQC directive syntax changed from assignment to directive-call form, but image strings are identical. No container image or dependency version was changed. |
| Resources and executor settings | Unchanged | CPUs, memory, time limits, executor settings, queues, and profiles were not edited. |
| Output paths and publishing | Unchanged | `publishDir` directives and output path patterns were not changed. |
| Sample parsing and filtering | Unchanged | CSV parsing, BAM/FASTQ channel construction, filters, grouping, combining, and branching logic were not intentionally changed. |

## Validation Summary

- `mamba run -n swgs_env nextflow lint .`: passed with warnings only.
- `mamba run -n swgs_env nextflow run main.nf -params-file config/params.yaml -profile singularity -preview`: passed; no tasks executed.
- `mamba run -n swgs_env env NXF_SYNTAX_PARSER=v1 nextflow run main.nf -params-file config/params.yaml -profile singularity -preview`: passed; no tasks executed.

## Human Review Items

- Confirm that `params.qd_bwgavgbed` is the intended bigWigAverageOverBed parameter for `QDNA_BINS`.
- Confirm WisecondorX ACN should consume the three outputs of `CN_WISECONDORX_RCN` in declared order: relative CN directory, bin size, BAM type.
- Decide separately whether to modernize deprecated `Channel` factory access and lint warnings; these were not changed because they are not strict-parser blockers.
