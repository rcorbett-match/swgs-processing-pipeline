# Nextflow 26 Strict Syntax Audit

## Baseline

- Repository: `/home/rcorbett/swgs-processing-pipeline`
- Baseline commit: `1681fc8`
- Baseline branch: `master`
- Baseline worktree status: clean
- Requested runtime: Nextflow `26.04.6`
- Validation environment: `mamba run -n swgs_env`
- Nextflow version in `swgs_env`: `26.04.6` build `12646`
- Java version in `swgs_env`: OpenJDK `25.0.2`
- Local default Nextflow outside `swgs_env`: `24.10.3`
- Local default Java outside `swgs_env`: OpenJDK `11.0.1`

Branch creation was attempted with `git switch -c nextflow-26-strict-syntax` before edits. The first attempt failed because `.git/refs` is read-only in the normal execution sandbox; the branch was then created successfully with escalated Git permissions.

## Inventory

Repository-wide Nextflow/config/parameter/test inventory:

- `main.nf`
- `workflows/swgs_qc.nf`
- `subworkflows/local/alignment_processing.nf`
- `subworkflows/local/cn_ichor.nf`
- `subworkflows/local/cn_qdnaseq.nf`
- `subworkflows/local/cn_wisecondorx.nf`
- `subworkflows/local/collect_bams.nf`
- `subworkflows/local/collect_fastqs.nf`
- `subworkflows/local/post_alignment_qc.nf`
- `subworkflows/local/pre_alignment_qc.nf`
- `modules/local/alignment/main.nf`
- `modules/local/azure_download/main.nf`
- `modules/local/cn_ichor/main.nf`
- `modules/local/cn_qdnaseq_acn/main.nf`
- `modules/local/cn_qdnaseq_rcn/main.nf`
- `modules/local/cn_wisecondorx_acn/main.nf`
- `modules/local/cn_wisecondorx_rcn/main.nf`
- `modules/local/fastqc_postalign/main.nf`
- `modules/local/fastqc_prealign/main.nf`
- `modules/local/ichor_combine/main.nf`
- `modules/local/ichor_process/main.nf`
- `modules/local/mark_duplicates/main.nf`
- `modules/local/multiqc_postalign/main.nf`
- `modules/local/multiqc_prealign/main.nf`
- `modules/local/qdna_bin_annot/main.nf`
- `modules/local/samtools_cov_flagstat/main.nf`
- `modules/local/samtools_sort_index/main.nf`
- `modules/local/trimmomatic/main.nf`
- `modules/local/wcx_ref_create/main.nf`
- `old_scripts/swgs-workflow-se-bwa-mem.nf`
- `old_scripts/swgs-workflow-se-test.nf`
- `nextflow.config`
- `old_scripts/nextflow_bwa-mem.config`
- `config/params.yaml`

No test directory was found.

## Baseline Diagnostics

Commands:

```bash
mamba run -n swgs_env nextflow -version
mamba run -n swgs_env java -version
mamba run -n swgs_env nextflow lint .
NXF_SYNTAX_PARSER=v1 nextflow run main.nf -params-file config/params.yaml -profile singularity -stub-run
```

Baseline `nextflow lint .` under Nextflow `26.04.6` failed with 73 errors in 9 files and 39 warnings in 10 files. Parser and compilation errors were grouped as follows:

- Undefined workflow call: `main.nf` included `SWGS_QC_PIPELINE` but called `WGS_QC_PIPELINE`.
- Top-level executable statements mixed with declarations: `workflows/swgs_qc.nf`, `subworkflows/local/post_alignment_qc.nf`, `old_scripts/swgs-workflow-se-test.nf`, and `old_scripts/swgs-workflow-se-bwa-mem.nf`.
- Invalid process directive syntax: `container = 'docker://staphb/fastqc:latest'` in both FastQC modules.
- Undefined interpolation: `modules/local/qdna_bin_annot/main.nf` referenced `qd_bwavgbed`, while `config/params.yaml` and `README.md` define `qd_bwgavgbed`.
- Incorrect process call arity: `subworkflows/local/cn_wisecondorx.nf` called `CN_WISECONDORX_RCN` and `CN_WISECONDORX_ACN` with too few arguments.

The baseline legacy-parser run under the local default `24.10.3` installation reached workflow execution and failed before tasks launched:

```text
Missing process or function WGS_QC_PIPELINE() -- Did you mean 'SWGS_QC_PIPELINE' instead?
```

## Changed Files

### `main.nf`

Original:

```nextflow
WGS_QC_PIPELINE()
```

Replacement:

```nextflow
SWGS_QC_PIPELINE()
```

Reason: the included workflow is named `SWGS_QC_PIPELINE`. This fixes an unambiguous name mismatch that prevented execution. It does not rename the workflow or change invocation order.

### `workflows/swgs_qc.nf`

Original construct: parameter assertions and parameter-derived assignments were executable top-level statements before `include` declarations and before `workflow SWGS_QC_PIPELINE`.

Replacement: the same assertion calls and assignments were moved, in the same order, to the beginning of `workflow SWGS_QC_PIPELINE`.

Reason: Nextflow 26 strict syntax rejects executable top-level statements mixed with script declarations. The expressions, parameter names, defaults, channel factory call, and branch order were preserved.

### `subworkflows/local/post_alignment_qc.nf`

Original:

```nextflow
multiqc_config = file("$projectDir/config/multiqc_config.yaml")
```

Replacement: the same assignment was moved into the `main:` block before `Channel.fromPath(multiqc_config, checkIfExists: true)`.

Reason: strict syntax rejects executable top-level assignment mixed with declarations. The file path and downstream channel construction were preserved.

### `modules/local/fastqc_prealign/main.nf`

Original:

```nextflow
container = 'docker://staphb/fastqc:latest'
```

Replacement:

```nextflow
container 'docker://staphb/fastqc:latest'
```

Reason: strict syntax rejects the assignment form as an invalid process directive. The container value is unchanged.

### `modules/local/fastqc_postalign/main.nf`

Original:

```nextflow
container = 'docker://staphb/fastqc:latest'
```

Replacement:

```nextflow
container 'docker://staphb/fastqc:latest'
```

Reason: strict syntax rejects the assignment form as an invalid process directive. The container value is unchanged.

### `modules/local/qdna_bin_annot/main.nf`

Original:

```nextflow
\$WORK_DIR/${qd_bwavgbed}
```

Replacement:

```nextflow
\$WORK_DIR/${params.qd_bwgavgbed}
```

Reason: `qd_bwavgbed` is undefined and prevents strict compilation. The configured parameter is `qd_bwgavgbed` in both `config/params.yaml` and `README.md`. This is a pre-existing name bug fixed because it blocks parsing/compilation.

### `subworkflows/local/cn_wisecondorx.nf`

Original:

```nextflow
CN_WISECONDORX_RCN(wx_input_ch)
CN_WISECONDORX_ACN(CN_WISECONDORX_RCN.out)
```

Replacement:

```nextflow
CN_WISECONDORX_RCN(wx_input_ch, params.output_directory)
CN_WISECONDORX_ACN(CN_WISECONDORX_RCN.out[0], CN_WISECONDORX_RCN.out[1], CN_WISECONDORX_RCN.out[2])
```

Original in the non-reference-creation branch:

```nextflow
CN_WISECONDORX_ACN(CN_WISECONDORX_RCN.out)
```

Replacement:

```nextflow
CN_WISECONDORX_ACN(CN_WISECONDORX_RCN.out[0], CN_WISECONDORX_RCN.out[1], CN_WISECONDORX_RCN.out[2])
```

Reason: strict lint validates process call arity. `CN_WISECONDORX_RCN` takes two inputs and `CN_WISECONDORX_ACN` takes three inputs. The replacement uses the existing output order from `CN_WISECONDORX_RCN`: relative CN directory, bin size, and BAM type. No channel operator order was changed.

### `old_scripts/swgs-workflow-se-test.nf`

Original construct: `log.info` and `reads_ch = Channel.fromPath(...).map(...)` were executable top-level statements.

Replacement: the same `log.info` and `reads_ch` construction were moved into the existing `workflow {}` block before `FASTQC(reads_ch)`.

Reason: strict syntax rejects executable top-level statements mixed with process declarations. Channel construction and workflow call order were preserved inside the old script.

### `old_scripts/swgs-workflow-se-bwa-mem.nf`

Original construct: `log.info` was an executable top-level statement.

Replacement: the same `log.info` block was moved into the existing `workflow {}` block before channel construction.

Reason: strict syntax rejects executable top-level statements mixed with process declarations. Existing workflow order was preserved.

## Suspected Pre-Existing Bugs or Non-Blocking Warnings Not Changed

- `modules/local/azure_download/main.nf` uses `projectDir` inside a process script. Nextflow 26 reports this as discouraged, not an error. It was left unchanged.
- Several files use `Channel` instead of `channel`. Nextflow 26 reports this as deprecated, not an error. These were left unchanged to avoid style-only changes.
- Several workflow `take:` variables are unused. These warnings were left unchanged because removing or renaming them would change workflow signatures.
- `subworkflows/local/collect_bams.nf` reports variables as declared but not used even though they are emitted. This was left unchanged.
- `subworkflows/local/collect_fastqs.nf` has a single unnamed emit warning. It was left unchanged.
- The `qd_bwavgbed` mismatch was a suspected pre-existing bug and was fixed only because strict lint treats it as an undefined symbol.
- The WisecondorX call arity mismatch was a suspected pre-existing bug and was fixed only because strict lint treats it as an error.

## Validation Results

Strict lint:

```bash
mamba run -n swgs_env nextflow lint .
```

Result: exit code 0. Repository-wide lint completed with 21 warnings and no errors.

Strict/default preview compilation:

```bash
mamba run -n swgs_env nextflow run main.nf -params-file config/params.yaml -profile singularity -preview
```

Result: exit code 0. The workflow compiled and printed the expected banner and parameter validation output. No tasks were executed because `-preview` was used.

Legacy parser preview comparison:

```bash
mamba run -n swgs_env env NXF_SYNTAX_PARSER=v1 nextflow run main.nf -params-file config/params.yaml -profile singularity -preview
```

Result: exit code 0. The same banner and parameter validation output were produced. No tasks were executed because `-preview` was used.

## Process Body and Behavior Preservation

All process input blocks, output blocks, resource directives, publish locations, labels, command tools, and command arguments were left unchanged except where listed above:

- FastQC container directive syntax changed, but the image string did not.
- QDNAseq bin annotation command interpolation changed from an undefined variable to the documented parameter name required for compilation.

No full production workflow was run. No production outputs were intentionally written or overwritten.

## Remaining Risks

- `-preview` validates compilation and workflow construction but does not prove scientific behavioral equivalence.
- No real or stub tasks were launched, so container availability and command runtime behavior were not tested.
- Existing lint warnings remain intentionally unchanged.
- Old scripts now parse under repository-wide lint, but they were not independently executed.
