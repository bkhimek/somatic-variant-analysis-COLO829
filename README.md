# COLO829 Tumour-Normal Somatic Genomics Pipeline

**Status:** Phase 0 (pre-implementation) — see `docs/PHASE0_FINDINGS.md` for the full research trail and `docs/data_sources.md` for the living data-source/version ledger.

Tumour-normal somatic variant calling (GATK Mutect2) benchmarked against a COLO829/COLO829BL truth set, with copy-number (CNVkit), mutational signature (SigProfilerMatrixGenerator + SigProfilerAssignment), and oncogenicity/actionability interpretation layers. Full design rationale: see the project plan (kept alongside this repo, not committed here — v1.1 as of 2026-08-27).

## Phase 0 outcome (2026-08-29)

Two things the plan asked to verify turned out to need correcting, not just confirming:

1. **Truth set citation was wrong.** The plan named SEQC2/Fang et al. 2021 — that paper benchmarks HCC1395, not COLO829. Real options and the decision needed: `docs/PHASE0_FINDINGS.md` §1.
2. **`SigProfilerSingleBase` doesn't exist on PyPI.** Correct packages: `SigProfilerMatrixGenerator` + `SigProfilerAssignment`. Detail: `docs/PHASE0_FINDINGS.md` §2.

Everything else on the Phase 0 checklist (FASTQ accessions, GATK hg38 resources, container images, COSMIC registration process, dev/full Nextflow profiles, run manifest schema) is resolved — see `docs/PHASE0_FINDINGS.md` for the full breakdown and what still needs your direct action (registration, actual downloads, confirming exact numbers from your own network access).

## Repo layout

```
COLO829-somatic-pipeline/
├── main.nf                  (Phase 1+)
├── nextflow.config           dev/full profiles + docker/singularity profiles (this batch)
├── modules/                  (Phase 1+, one process per file)
├── workflows/                (Phase 1+)
├── bin/                      (Phase 1+)
├── conf/resources.config     configurable QC thresholds (this batch)
├── data/gene_lists/          melanoma driver gene seed list (this batch)
├── docs/
│   ├── PHASE0_FINDINGS.md    Phase 0 research report (this batch)
│   ├── data_sources.md       living data-source/version/licensing ledger (this batch)
│   ├── run_manifest.schema.json + run_manifest.example.json  (this batch)
│   ├── somatic_interpretation.md   (Phase 5+)
│   └── benchmarking_results.md     (Phase 3+)
└── .gitignore                 (this batch)
```

## Syncing into the real repo

Built in a Cowork sandbox (no direct WSL access this session). To land it in `~/projects/somatic-variant-analysis-COLO829/`:

```bash
cd ~/projects/somatic-variant-analysis-COLO829
unzip -o /path/to/phase0_batch.zip -d .
git status
git diff --stat
git add -A
git commit -m "Phase 0: pre-implementation research, corrected truth-set citation and SigProfiler package name, dev/full Nextflow profiles, run manifest schema"
git push
```
