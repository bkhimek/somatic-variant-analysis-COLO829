# COLO829 Tumour-Normal Somatic Genomics Pipeline

**Status:** Phase 1 (QC + alignment) and Phase 2 (contamination + Mutect2 somatic calling) both signed off; Phase 3 (benchmarking against the NYGC truth set) built and wired in, not yet run — see `docs/PHASE1_NOTES.md`, `docs/PHASE2_NOTES.md`, and `docs/PHASE3_NOTES.md` for what each phase actually proved (and explicitly didn't), `docs/PHASE0_FINDINGS.md` for the Phase 0 research trail, and `docs/data_sources.md` for the living data-source/version ledger.

Tumour-normal somatic variant calling (GATK Mutect2) benchmarked against a COLO829/COLO829BL truth set, with copy-number (CNVkit), mutational signature (SigProfilerMatrixGenerator + SigProfilerAssignment), and oncogenicity/actionability interpretation layers. Full design rationale: see the project plan (kept alongside this repo, not committed here — v1.1 as of 2026-08-27).

## Phase 0 outcome (2026-08-29)

Two things the plan asked to verify turned out to need correcting, not just confirming:

1. **Truth set citation was wrong.** The plan named SEQC2/Fang et al. 2021 — that paper benchmarks HCC1395, not COLO829. **Decided:** NYGC's open COLO-829 (HiSeqX) somatic VCF is the primary SNV/indel truth set (Craig et al. 2016 would be the more authoritative multi-platform consensus, but is confirmed EGA/dbGaP controlled-access only with no open alternative — declined to preserve the project's controlled-access-free scope). Full rationale and caveats: `docs/PHASE0_FINDINGS.md` §1, `docs/data_sources.md` §3.
2. **`SigProfilerSingleBase` doesn't exist on PyPI.** Correct packages: `SigProfilerMatrixGenerator` + `SigProfilerAssignment`. Detail: `docs/PHASE0_FINDINGS.md` §2.

Everything else on the Phase 0 checklist (FASTQ accessions, truth set — decided, downloaded, extracted — GATK hg38 resources, container images, COSMIC registration process, dev/full Nextflow profiles, run manifest schema) is resolved — see `docs/PHASE0_FINDINGS.md` for the full breakdown.

## Phase 1 status (2026-08-30) — SIGNED OFF

Modules 1–2 (QC, Alignment) ran successfully end-to-end against the real full genome, with ~99%/98.8% mapped fractions confirming the pipeline is correct. Getting there took five real bugs found only by actually executing the code (three Nextflow 26.x compatibility breaks, one resource-sizing correction, one cosmetic handler issue) — full story in `docs/PHASE1_NOTES.md`.

## Phase 2 status (2026-08-30) — SIGNED OFF

Modules 3–4 (Contamination estimation, Mutect2 somatic calling) ran successfully end-to-end on the `dev` profile (driver-gene-BED-restricted), confirming `GetPileupSummaries` → `CalculateContamination` → `Mutect2` → `LearnReadOrientationModel` → `FilterMutectCalls` are all wired correctly against the real reference and real GATK resource bundles. Getting there took two real bugs found only by execution (`GetPileupSummaries`'s mandatory `-L`, a GATK/JVM-in-Docker heap-sizing issue) plus one real finding that isn't a bug: research into GATK's own guidance, Broad's production WDL, and nf-core/sarek confirmed unsharded genome-wide Mutect2 execution is not how this tool is meant to run anywhere, at any memory size available on this machine. Real interval-scatter/gather architecture (and likely another one-off cloud compute burst) are explicitly deferred to when real full-depth data makes them actually necessary — full story in `docs/PHASE2_NOTES.md`.

## Phase 3 status (2026-08-30) — built, not yet run

Module 5 (benchmarking, via hap.py against the NYGC COLO829 truth set) is implemented in `modules/benchmarking.nf`, extending the same `workflows/somatic.nf`. Since the only Mutect2 output that exists so far has zero variant calls (Phase 2's `dev`-profile result), this first run is only expected to prove the DAG/container/hap.py command line are wired correctly — not produce meaningful precision/recall numbers, which stay deferred until real full-genome Mutect2 execution happens. **Read `docs/PHASE3_NOTES.md` before your first Phase 3 run** — it covers the new required param, the design choices behind the hap.py invocation, and a real (if low-risk) contig-naming check worth doing if it errors out.

## Repo layout

```
COLO829-somatic-pipeline/
├── main.nf                   Entry point (Phase 1 + Phase 2 + Phase 3)
├── nextflow.config            dev/full profiles + docker/singularity profiles
├── assets/NO_FILE             optional-input sentinel (Phase 2)
├── modules/
│   ├── fastqc.nf              Module 1: FASTQC + MULTIQC
│   ├── alignment.nf           Module 2: BWA-MEM2 -> sort -> MarkDuplicates
│   ├── reference_prep.nf      Phase 2: .fai/.dict/VCF-index auto-prep
│   ├── contamination.nf       Module 3: GetPileupSummaries + CalculateContamination
│   ├── mutect2.nf             Module 4: Mutect2 + LearnReadOrientationModel + FilterMutectCalls
│   └── benchmarking.nf        Module 5: PrepareTruthVcf + hap.py benchmarking
├── workflows/somatic.nf        wires Modules 1-5 together
├── bin/                       (Phase 5+)
├── conf/resources.config      QC + contamination thresholds (includeConfig'd from Phase 2)
├── data/gene_lists/           melanoma driver gene seed list + dev_intervals.bed (real GRCh38 coords)
├── docs/
│   ├── PHASE0_FINDINGS.md     Phase 0 research report
│   ├── PHASE1_NOTES.md        Phase 1 implementation notes (signed off)
│   ├── PHASE2_NOTES.md        Phase 2 implementation notes (signed off)
│   ├── PHASE3_NOTES.md        Phase 3 implementation notes -- read before running
│   ├── data_sources.md        living data-source/version/licensing ledger
│   ├── run_manifest.schema.json + run_manifest.example.json
│   ├── somatic_interpretation.md   (Phase 5+)
│   └── benchmarking_results.md     (written once real full-genome results exist)
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
