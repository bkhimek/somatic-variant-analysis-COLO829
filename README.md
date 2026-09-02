# COLO829 Tumour-Normal Somatic Genomics Pipeline

**New to this repo? See [`PORTFOLIO_SUMMARY.md`](PORTFOLIO_SUMMARY.md) for a top-level narrative of what this project demonstrates and what it actually proved — this README is the ongoing technical status ledger, updated phase by phase as work happens.**

**Status:** Phase 1 (QC + alignment), Phase 2 (contamination + Mutect2 somatic calling), Phase 3 (benchmarking against the NYGC truth set), and Phase 4 (CNVkit copy-number calling) all signed off. Module 4 (Mutect2) has since been revisited (2026-09-01) to close out the "unsharded genome-wide Mutect2 doesn't fit this machine" limitation Phase 2 deferred — real interval scatter/gather, designed against GATK's own docs, Broad's production WDL, and nf-core/sarek's source, built and confirmed via a real 2026-09-02 run against real data. Module 8 (oncogenicity/actionability interpretation, CIViC-based) is **SIGNED OFF** as of 2026-09-02 — confirmed via a real pipeline run against the real data, matching an independent standalone test exactly. See `docs/ONCOGENICITY_NOTES.md`. See `docs/PHASE1_NOTES.md`, `docs/PHASE2_NOTES.md`, `docs/PHASE3_NOTES.md`, `docs/PHASE4_NOTES.md`, and `docs/MUTECT2_SCATTERGATHER_NOTES.md` for what's actually been proved (and explicitly hasn't), `docs/PHASE0_FINDINGS.md` for the Phase 0 research trail, `docs/REAL_DATA_RESULTS.md` for the first real-data findings, `docs/ONCOGENICITY_NOTES.md` for Module 8, and `docs/data_sources.md` for the living data-source/version ledger.

Tumour-normal somatic variant calling (GATK Mutect2) benchmarked against a COLO829/COLO829BL truth set, with copy-number (CNVkit), mutational signature (SigProfilerMatrixGenerator + SigProfilerAssignment), and oncogenicity/actionability interpretation layers. Full design rationale: see the project plan (kept alongside this repo, not committed here — v1.1 as of 2026-08-27).

## Phase 0 outcome (2026-08-29)

Two things the plan asked to verify turned out to need correcting, not just confirming:

1. **Truth set citation was wrong.** The plan named SEQC2/Fang et al. 2021 — that paper benchmarks HCC1395, not COLO829. **Decided:** NYGC's open COLO-829 (HiSeqX) somatic VCF is the primary SNV/indel truth set (Craig et al. 2016 would be the more authoritative multi-platform consensus, but is confirmed EGA/dbGaP controlled-access only with no open alternative — declined to preserve the project's controlled-access-free scope). Full rationale and caveats: `docs/PHASE0_FINDINGS.md` §1, `docs/data_sources.md` §3.
2. **`SigProfilerSingleBase` doesn't exist on PyPI.** Correct packages: `SigProfilerMatrixGenerator` + `SigProfilerAssignment`. Detail: `docs/PHASE0_FINDINGS.md` §2.

Everything else on the Phase 0 checklist (FASTQ accessions, truth set — decided, downloaded, extracted — GATK hg38 resources, container images, COSMIC registration process, dev/full Nextflow profiles, run manifest schema) is resolved — see `docs/PHASE0_FINDINGS.md` for the full breakdown.

## Phase 1 status (2026-08-30) — SIGNED OFF

Modules 1–2 (QC, Alignment) ran successfully end-to-end using the full GRCh38 reference genome as the alignment target, with ~99%/98.8% mapped fractions confirming the pipeline is correct. **Clarified 2026-09-01** (see `docs/MUTECT2_SCATTERGATHER_NOTES.md`): this line previously read "against the real full genome" in a way that could be misread as real full sequencing depth — `docs/PHASE1_NOTES.md`'s own record is precise ("Same 10,000-read-pair subsample ... this time aligned against the real full genome" i.e. the reference, not the input depth), but this summary line wasn't, and that ambiguity led to a real wrong assumption later in the project. No real full-depth BAM has ever actually existed — every dedup BAM produced through 2026-09-01 is this same 10,000-read-pair subsample. Getting Modules 1-2 working took five real bugs found only by actually executing the code (three Nextflow 26.x compatibility breaks, one resource-sizing correction, one cosmetic handler issue) — full story in `docs/PHASE1_NOTES.md`.

## Phase 2 status (2026-08-30) — SIGNED OFF

Modules 3–4 (Contamination estimation, Mutect2 somatic calling) ran successfully end-to-end on the `dev` profile (driver-gene-BED-restricted), confirming `GetPileupSummaries` → `CalculateContamination` → `Mutect2` → `LearnReadOrientationModel` → `FilterMutectCalls` are all wired correctly against the real reference and real GATK resource bundles. Getting there took two real bugs found only by execution (`GetPileupSummaries`'s mandatory `-L`, a GATK/JVM-in-Docker heap-sizing issue) plus one real finding that isn't a bug: research into GATK's own guidance, Broad's production WDL, and nf-core/sarek confirmed unsharded genome-wide Mutect2 execution is not how this tool is meant to run anywhere, at any memory size available on this machine. Real interval-scatter/gather architecture (and likely another one-off cloud compute burst) are explicitly deferred to when real full-depth data makes them actually necessary — full story in `docs/PHASE2_NOTES.md`.

## Phase 3 status (2026-08-30) — SIGNED OFF

Module 5 (benchmarking against the NYGC COLO829 truth set) is implemented in `modules/benchmarking.nf`, extending the same `workflows/somatic.nf`, and completed its first successful end-to-end run against real data. Uses **som.py**, not hap.py, despite both shipping in the same container — hap.py expects a GT (genotype) field to judge whether a record is a real call, and neither the NYGC truth set nor Mutect2's own output carry one (both use AD/DP/AF-style somatic annotation instead), which silently zeroed out every real variant and made hap.py hard-fail; som.py is Illumina's own somatic-specific comparison tool built for exactly this. The truth set itself also had to change: the originally-chosen COLO-829 HiSeqX VCF turned out to be corrupted at the source in NYGC's own published archive (confirmed via byte-level forensics, with no fix available anywhere) — now using the COLO-829 NovaSeq VCF from the same archive instead, which reopens a platform-mismatch caveat (our FASTQs are HiSeq X Ten). Since the only Mutect2 output that exists so far has zero variant calls (Phase 2's `dev`-profile result), the successful run proved the DAG/container/som.py command line are wired correctly (`som.stats.csv` came back with the exact predicted 0 TP/0 FP/all-FN result) — it did not, and was never meant to, produce meaningful precision/recall numbers, which stay deferred until real full-genome Mutect2 execution happens. Signed off on that basis — full story, including what's explicitly deferred, in `docs/PHASE3_NOTES.md`'s "Phase 3 sign-off" section.

## Phase 4 status (2026-08-31) — SIGNED OFF

Module 6 (CNVkit whole-genome tumour/normal copy-number calling) is implemented in `modules/cnvkit.nf`, extending the same `workflows/somatic.nf`, and completed a clean end-to-end run against real data. Designed directly from CNVkit's own documentation (the project plan's exact CNVkit wording wasn't retrieved for this build — see `docs/PHASE4_NOTES.md` for that scoping decision) rather than guessed: whole-genome mode (`-m wgs`, no target/antitarget BED needed — CNVkit computes accessible regions on the fly), `-y`/`--male-reference` included because the COLO829 donor's sex was independently verified (male, per ATCC's own CRL-1974 product page) rather than left at CNVkit's female-reference default, and `--drop-low-coverage` included per CNVkit's own recommendation for tumour samples. `--annotate`/`--scatter`/`--diagram` are deliberately deferred for this first cut. Getting to a clean run took three real findings, each root-caused by execution rather than guessed: CNVkit's auto-computed ~53kb genome-wide bins were far too small for the dev-profile subsample's ~0.001x effective coverage (100% of bins failed CNVkit's own reference-quality filter), fixed with a `dev`-profile-only 10Mb bin-size override; `batch` wrote three `.cns`-suffixed files rather than the one this module first assumed, breaking the hand-off to `call`, fixed with exact output filenames and separated `publishDir`s; and Nextflow's `-resume` turned out not to re-derive a cached task's outputs from a changed `output:` block, which delayed confirming the second fix until that task's work directory was deleted to force a real re-execution. The resulting `.call.cns` is exactly the "mechanically correct, no real signal" result expected from this near-empty subsample (5 single-probe segments, degenerate confidence intervals, depths of a few reads per 10Mb) — real copy-number findings stay deferred to real full-genome data, same basis as Phase 2/3's sign-offs. Full story in `docs/PHASE4_NOTES.md`.

## Real-data run (2026-09-02) — SIGNED OFF, real driver mutations found and cross-checked

The gene-panel real-data extraction (`bin/extract_real_gene_panel.sh`, `docs/data_sources.md` §1) succeeded with sane, proportional read counts, and the resulting real FASTQs were run through the full pipeline end-to-end for the first time (26/26 processes succeeded). Mutect2 reported 8 `PASS` calls; cross-checked against [Cellosaurus's independently curated COLO829 genotype (CVCL_1137)](https://www.cellosaurus.org/CVCL_1137), 3 of them are exact matches to COLO829's documented driver mutations — **BRAF V600E** (chr7:140753336 A>T, the canonical hg38 hotspot, confirmed independently via OncoKB/GeneBe), a **CDKN2A frameshift deletion** (chr9:21971154, matching the documented `c.203_204delCG`), and a **TERT promoter dinucleotide mutation** (chr5:1295113 GG>AA, the reverse-complement of the documented `CC>TT` promoter change) — and `som.py`'s benchmarking against the independent NYGC truth set confirms **all 8 calls as true positives, zero false positives**. A documented, known-scope gap (PTEN's large 142bp deletion, not recoverable by an SNV/indel caller) and a documented reason CNVkit's genome-wide output isn't interpretable here (real coverage exists only over the gene panel) are both disclosed rather than glossed over. Full detail, every coordinate, and every source checked: `docs/REAL_DATA_RESULTS.md`.

## Mutect2 interval scatter/gather (2026-09-01) — dev-profile DAG confirmed, full-profile not yet run

Phase 2 deferred real full-genome Mutect2 execution after an unsharded, genome-wide invocation didn't fit this machine at any memory size tried, pending real full-depth data making it necessary. `modules/mutect2.nf`'s `MUTECT2` process is now scattered across `SplitIntervals`-generated interval shards (`ScatterIntervalsByNs` first on the `full` profile, to keep shard boundaries out of assembly gaps) instead of running once against the whole genome, with `MergeVcfs`/`MergeMutectStats`/a multi-input `LearnReadOrientationModel` call gathering the per-shard outputs back into exactly the same one-per-sample shape `FilterMutectCalls` always expected. Every choice here — which merge tool, whether f1r2 tarballs get pre-merged, subdivision mode — was checked directly against GATK's own documentation, Broad's production `mutect2.wdl`, and nf-core/sarek's actual Nextflow source, not guessed. A `dev`-profile smoke test (2 shards) confirmed the DAG end-to-end: `SplitIntervals` → 2x `MUTECT2` → the three merge/gather steps → `FilterMutectCalls` → `SOMPY_BENCHMARK` all completed, with output identical (0 TP/0 FP/all-FN) to every prior unsharded run against the same sparse dev data. The scatter count itself (`params.mutect2_scatter_count`: 2 on `dev`, 20 on `full`) is an explicitly-flagged starting guess for `full`.

**Correction, 2026-09-01:** this section originally said real full-genome BAMs already existed from Phase 1, justifying jumping straight to a real `full`-profile run. That was wrong — checking actual read counts found every dedup BAM this project has ever produced is the same 10,000-read-pair subsample; no real full-depth data has ever existed here (the README's Phase 1 line was ambiguous in a way that caused this — see that section above). **Real path forward, found instead:** ENA's PRJEB27698 submission (`docs/data_sources.md` §1) includes pre-aligned, deduplicated BAMs for both real COLO829/COLO829BL sequencing runs, each with a `.bai` index — but aligned to GRCh37, not this pipeline's GRCh38. Rather than a full ~174GB download and re-alignment, the plan is to extract real reads over the melanoma gene panel directly from those remote BAMs (`samtools view -L`, using GRCh37 coordinates for the same genes in the new `data/gene_lists/dev_intervals_grch37.bed`), convert back to FASTQ, and re-align through this pipeline's own Module 1-2 to GRCh38 — real ~37X/~98X depth in the panel region, without the full-genome download/compute commitment. See `docs/data_sources.md` §1 for the extraction commands and `docs/MUTECT2_SCATTERGATHER_NOTES.md` for the full story. **Not yet run** — extraction commands are documented but haven't been executed.

## Repo layout

```
COLO829-somatic-pipeline/
├── main.nf                   Entry point (Phase 1 + Phase 2 + Phase 3 + Phase 4)
├── nextflow.config            dev/full profiles + docker/singularity profiles
├── assets/NO_FILE             optional-input sentinel (Phase 2)
├── modules/
│   ├── fastqc.nf              Module 1: FASTQC + MULTIQC
│   ├── alignment.nf           Module 2: BWA-MEM2 -> sort -> MarkDuplicates
│   ├── reference_prep.nf      Phase 2: .fai/.dict/VCF-index auto-prep
│   ├── contamination.nf       Module 3: GetPileupSummaries + CalculateContamination
│   ├── mutect2.nf             Module 4: Mutect2 + LearnReadOrientationModel + FilterMutectCalls
│   ├── benchmarking.nf        Module 5: PrepareTruthVcf + som.py benchmarking
│   ├── cnvkit.nf              Module 6: CNVkit whole-genome tumour/normal copy-number calling
│   └── oncogenicity.nf        Module 8: PASS-filter -> SnpEff hg38 annotation -> CIViC lookup (docs/ONCOGENICITY_NOTES.md)
├── workflows/somatic.nf        wires Modules 1-8 together
├── bin/
│   ├── extract_real_gene_panel.sh  robust per-gene retry extraction of real reads from ENA's
│   │                                pre-aligned GRCh37 BAMs (added 2026-09-02, see docs/data_sources.md §1)
│   └── civic_annotate.py      Module 8: matches SnpEff HGVS.p output to CIViC evidence (docs/ONCOGENICITY_NOTES.md)
├── conf/resources.config      QC + contamination thresholds (includeConfig'd from Phase 2)
├── data/gene_lists/           melanoma driver gene seed list + dev_intervals.bed (GRCh38 coords)
│                              + dev_intervals_grch37.bed (same genes, GRCh37 -- for extracting real reads
│                              from ENA's pre-aligned GRCh37 BAMs, see docs/data_sources.md §1)
├── docs/
│   ├── PHASE0_FINDINGS.md     Phase 0 research report
│   ├── PHASE1_NOTES.md        Phase 1 implementation notes (signed off)
│   ├── PHASE2_NOTES.md        Phase 2 implementation notes (signed off)
│   ├── PHASE3_NOTES.md        Phase 3 implementation notes (signed off)
│   ├── PHASE4_NOTES.md        Phase 4 implementation notes (signed off)
│   ├── MUTECT2_SCATTERGATHER_NOTES.md  Mutect2 interval scatter/gather -- built, not yet run
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
