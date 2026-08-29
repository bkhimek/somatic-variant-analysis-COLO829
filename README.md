# COLO829 Tumour-Normal Somatic Genomics Pipeline

**Status:** Phase 1 (QC + alignment) written, not yet run on real data — see `docs/PHASE1_NOTES.md` for what to check before your first `nextflow run`, `docs/PHASE0_FINDINGS.md` for the Phase 0 research trail, and `docs/data_sources.md` for the living data-source/version ledger.

Tumour-normal somatic variant calling (GATK Mutect2) benchmarked against a COLO829/COLO829BL truth set, with copy-number (CNVkit), mutational signature (SigProfilerMatrixGenerator + SigProfilerAssignment), and oncogenicity/actionability interpretation layers. Full design rationale: see the project plan (kept alongside this repo, not committed here — v1.1 as of 2026-08-27).

## Phase 0 outcome (2026-08-29)

Two things the plan asked to verify turned out to need correcting, not just confirming:

1. **Truth set citation was wrong.** The plan named SEQC2/Fang et al. 2021 — that paper benchmarks HCC1395, not COLO829. **Decided:** NYGC's open COLO-829 (HiSeqX) somatic VCF is the primary SNV/indel truth set (Craig et al. 2016 would be the more authoritative multi-platform consensus, but is confirmed EGA/dbGaP controlled-access only with no open alternative — declined to preserve the project's controlled-access-free scope). Full rationale and caveats: `docs/PHASE0_FINDINGS.md` §1, `docs/data_sources.md` §3.
2. **`SigProfilerSingleBase` doesn't exist on PyPI.** Correct packages: `SigProfilerMatrixGenerator` + `SigProfilerAssignment`. Detail: `docs/PHASE0_FINDINGS.md` §2.

Everything else on the Phase 0 checklist (FASTQ accessions, truth set — decided, downloaded, extracted — GATK hg38 resources, container images, COSMIC registration process, dev/full Nextflow profiles, run manifest schema) is resolved — see `docs/PHASE0_FINDINGS.md` for the full breakdown.

## Phase 1 status (2026-08-29)

Modules 1–2 (QC, Alignment) are implemented in `modules/`, `workflows/somatic.nf`, and `main.nf`, following the plan's module spec. **Not yet run** — built without Nextflow/Docker/data access in this sandbox, so it's been checked by manual DSL2 review only (which did catch and fix two real channel-cardinality bugs) rather than an actual execution. **Read `docs/PHASE1_NOTES.md` before your first run** — it lists an unconfirmed container tag, a design note on why bwa-mem2 and samtools are separate processes, and an important limitation on what the `dev` profile does and doesn't speed up.

## Repo layout

```
COLO829-somatic-pipeline/
├── main.nf                   Phase 1 entry point (this batch)
├── nextflow.config            dev/full profiles + docker/singularity profiles
├── modules/
│   ├── fastqc.nf              Module 1: FASTQC + MULTIQC (this batch)
│   └── alignment.nf           Module 2: BWA-MEM2 -> sort -> MarkDuplicates (this batch)
├── workflows/somatic.nf        wires Modules 1-2 together (this batch)
├── bin/                       (Phase 5+)
├── conf/resources.config      configurable QC thresholds
├── data/gene_lists/           melanoma driver gene seed list
├── docs/
│   ├── PHASE0_FINDINGS.md     Phase 0 research report
│   ├── PHASE1_NOTES.md        Phase 1 implementation notes -- read before running (this batch)
│   ├── data_sources.md        living data-source/version/licensing ledger
│   ├── run_manifest.schema.json + run_manifest.example.json
│   ├── somatic_interpretation.md   (Phase 5+)
│   └── benchmarking_results.md     (Phase 3+)
└── .gitignore
```

## Syncing into the real repo

Built in a Cowork sandbox (no direct WSL access this session). The zip contains one wrapper folder (`COLO829-somatic-pipeline/`) that needs flattening into the repo root — don't unzip straight into the repo dir. To land an update batch in `~/projects/somatic-variant-analysis-COLO829/`:

```bash
ZIP="/mnt/c/Users/krist/OneDrive/Documents/Projects/Somatic_Variant_Analysis_Pipeline_COLO829_Malanoma/zip_files/<batch>.zip"
rm -rf /tmp/colo829_extract && mkdir -p /tmp/colo829_extract
unzip -o "$ZIP" -d /tmp/colo829_extract

cd ~/projects/somatic-variant-analysis-COLO829
cp -r /tmp/colo829_extract/COLO829-somatic-pipeline/. .
rm -rf /tmp/colo829_extract

git status
git diff --stat
git add -A
git commit -m "<describe what changed in this batch>"
git push
```

**Repo status as of 2026-08-29:** initialized, remote set to `https://github.com/bkhimek/somatic-variant-analysis-COLO829.git`, first commit pushed to `main`.
