# Phase 0 — Pre-Implementation Findings
**Project 8 — COLO829 Tumour-Normal Somatic Genomics Pipeline**
**Date:** 2026-08-29
**Status:** Research complete; two items need your decision before Phase 1 starts (flagged ⚠️ below)

This document reports what was resolved against the Phase 0 checklist in the v1.1 plan (§3 and §11). Full citations and download commands are in `docs/data_sources.md`. Items requiring your direct action (registration, actual downloads, running commands against your own network) are marked accordingly — none of that can be done from this sandbox.

---

## ⚠️ 1. Truth set — the plan's citation needs correcting

The plan names "SEQC2 / Fang et al., *Nature Biotechnology* 2021" as the COLO829 truth set, flagged "to be verified." It doesn't hold up:

**Fang et al. 2021 (the actual SEQC2 Somatic Mutation Working Group flagship paper) benchmarks HCC1395/HCC1395BL, not COLO829.** Its released call sets, at `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/`, are for HCC1395. COLO829 doesn't appear in that paper or on the SEQC2 consortium site (sites.google.com/view/seqc2) as one of the working group's reference pairs. This is exactly the kind of thing Phase 0 exists to catch before it propagates into `docs/benchmarking_results.md`.

The real options for a COLO829 SNV/indel truth set are less clean than the plan assumed:

| Source | Scope | Access | Genome build | Fit |
|---|---|---|---|---|
| **Craig et al. 2016**, *Sci Rep* (["A somatic reference standard for cancer genome sequencing"](https://www.nature.com/articles/srep24607)) | SNV + indel + CNV, 3-platform consensus (>35,000 point mutations, 446 indels) | Raw data on **EGA, controlled access** (study [EGAS00001001385](https://ega-archive.org/studies/EGAS00001001385), DAC [EGAC00001000408](https://www.omicsdi.org/dataset/ega/EGAC00001000408)) — the *call set itself* may be available as a paper supplementary file even though the reads require DAC approval; needs direct checking against the paper's supplementary materials | Not confirmed (2016-era, likely GRCh37 originally) | Widely cited as "the" COLO829 gold standard, but not simply "open access" as the plan assumed |
| **Valle-Inclan et al. 2022**, *Cell Genomics* (["A multi-platform reference for somatic structural variation detection"](https://www.sciencedirect.com/science/article/pii/S2666979X22000726)) | **Structural variants only** — no SNV/indel truth set | Fully open, ENA [PRJEB27698](https://www.ebi.ac.uk/ena/browser/view/PRJEB27698), truth VCFs on Zenodo ([4716169](https://zenodo.org/records/4716169), GRCh37 + GRCh38-lifted) | GRCh37 native, GRCh38-lifted VCF also provided | Same sample lineage as the open FASTQs below, but doesn't cover what Module 5 needs (SNV/indel) |
| **NYGC / Xiao et al. 2019**, *Sci Rep* (["Deep whole-genome sequencing of 3 cancer cell lines on 2 sequencing platforms"](https://www.nature.com/articles/s41598-019-55636-3)) | SNV/indel/CNV/SV somatic calls, HiSeqX + NovaSeq, up to 278X | **Processed somatic VCFs open, no login**, via the [NYGC companion site](https://www.nygenome.org/bioinformatics/3-cancer-cell-lines-on-2-sequencers/); raw FASTQ/BAM on **dbGaP phs001839 (controlled access)** | GRCh38 | Open SNV/indel calls, but from NYGC's own pipeline (not an independent multi-platform consensus like Craig 2016) and from different library preps than PRJEB27698 |

**None of these is a drop-in replacement that matches the plan's "no controlled-access issues" framing for a truth set.** The cleanest path given the "no controlled-access" constraint is: use the NYGC open somatic VCF as the SNV/indel truth set, explicitly document that it comes from different sequencing (HiSeqX/NovaSeq, NYGC library prep) than the ENA FASTQs you'll actually align — this is precisely the "cell lines evolve, truth set may not match preparation" caveat the plan already told us to document in `docs/benchmarking_results.md` (§11, Phase 0 truth-set bullet). Optionally cross-reference against Craig et al. 2016 supplementary calls if those turn out to be open (unverified — see TODO below).

**This is a decision for you, not something I should pick unilaterally** — it changes what "benchmarking success" means for the project's centrepiece deliverable. Options:
- (a) NYGC open VCF as primary truth set, caveat documented, Craig 2016 as a secondary sanity-check if its supplementary calls turn out to be open
- (b) Pursue EGA DAC access for Craig et al. 2016 (breaks the "no controlled-access issues" framing that motivated choosing COLO829 over TCGA in the first place)
- (c) Restrict truth-set validated scope to structural variants (Valle-Inclan, fully open) and treat SNV/indel benchmarking as best-effort against NYGC calls only

I'd lean (a), but this is exactly the kind of call the plan says should gate downstream work, so flagging rather than deciding for you.

---

## ⚠️ 2. Tool name doesn't exist: `SigProfilerSingleBase`

`pip install SigProfilerSingleBase` fails — **no such package on PyPI.** Tested directly in this session:
```
ERROR: Could not find a version that satisfies the requirement SigProfilerSingleBase
```
The Alexandrov Lab's real, current packages for this job:
- **`SigProfilerMatrixGenerator`** (PyPI, latest 1.3.6) — builds the 96-channel SBS trinucleotide matrix from a VCF
- **`SigProfilerAssignment`** (PyPI, latest 1.1.5) — decomposes/attributes an observed spectrum against a reference signature set (COSMIC SBS v3.4) — this is the actual tool for what Module 9 wants ("decomposed against COSMIC SBS reference set v3.4")
- `SigProfilerExtractor` exists too (latest 1.4.1) but does *de novo* extraction, which isn't what Module 9 asks for (it wants attribution against a known reference, not discovering new signatures)

**Action:** replace `SigProfilerSingleBase` everywhere in the plan/repo with `SigProfilerMatrixGenerator` + `SigProfilerAssignment` in `requirements.txt` and Module 9's spec. This is a factual correction, not a design change — the module's actual intent (decompose against COSMIC v3.4 reference) is unaffected.

---

## 3. FASTQ accessions — resolved

**ENA project [PRJEB27698](https://www.ebi.ac.uk/ena/browser/view/PRJEB27698)** ([Valle-Inclan et al. 2022](https://www.biorxiv.org/content/10.1101/2020.10.15.340497v1.full) data availability statement) is confirmed **fully open access, no approval required**, and includes Illumina HiSeq X Ten short-read WGS for both COLO829 (tumour) and COLO829BL (normal) — this is the same project that produced the open SV truth set above, so FASTQ and (SV) truth set share lineage, which is a plus.

**Caveat:** PRJEB27698 is a *multi-platform* project — 29 SRA experiments, ~4.09 Tbases / 1.55 TB total, spanning Illumina HiSeq, 10x Genomics linked-reads, Oxford Nanopore, PacBio, and BioNano optical mapping. **Do not bulk-download the whole project** — only the plain Illumina HiSeq X Ten WGS runs are relevant to this pipeline's BWA-MEM2/Mutect2 design.

**Outstanding (needs to run from your machine, not this sandbox):** I could not get past ENA's API rate limit from this sandbox to pull the exact run-accession-level table (which ERR runs are the plain Illumina WGS ones, vs. 10x/nanopore/PacBio/BioNano). Run this from WSL, where you have unrestricted network access:
```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB27698&result=read_run&fields=run_accession,sample_title,instrument_platform,library_strategy,fastq_ftp,base_count&format=tsv" \
  | awk -F'\t' '$3=="ILLUMINA"'
```
This gives you the run accessions, FASTQ FTP links, and base counts needed to fill in the compute/storage estimate (§5 below) with real numbers instead of the rule-of-thumb estimate I've given.

---

## 4. GATK hg38 reference resources — resolved

From the Broad's public GCS bucket (`gs://gatk-best-practices/somatic-hg38/`), confirmed via [GATK's Resource Bundle docs](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) and community threads:

| Resource | Path |
|---|---|
| Panel of Normals | `gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz` |
| gnomAD germline resource (Mutect2 `--germline-resource`) | `gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz` |
| Common biallelic sites (GetPileupSummaries, `--variant` for contamination estimation) | `gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz` |

Also downloadable without a GCP account via the public HTTPS mirror: `https://storage.googleapis.com/gatk-best-practices/somatic-hg38/<filename>`.

**Note on gnomAD version:** this bucket's `af-only-gnomad.hg38.vcf.gz` is the Broad's own curated Mutect2-ready extract, not necessarily the literal current gnomAD v4.1 release the plan cites in §8. Document the actual embedded gnomAD version once downloaded (check the VCF header) — this is exactly what `run_manifest.json`'s `gnomad_version` field (§5 below) is for.

---

## 5. Container images — one correction, one caveat

- **CNVkit:** `quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0` as specified in the plan **is a valid, real bioconda/biocontainers build** (confirmed 0.9.10-0 exists in the bioconda recipe index). Direct tag-list confirmation was blocked by quay.io's robots.txt from this sandbox — recommend a quick `docker pull quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0` on your machine before wiring it into `cnvkit.nf`, since biocontainers build-hash suffixes occasionally shift between bioconda rebuilds.
- **hap.py:** the plan doesn't pin a tag beyond "`pkrusche/hap.py:latest` (or equivalent — verify image availability)." Verified: **`pkrusche/hap.py:latest` was last pushed ~9 years ago** — it still pulls, but it's unmaintained. Two better-maintained alternatives, both confirmed to exist:
  - `quay.io/biocontainers/hap.py:0.3.15-0` (bioconda's current hap.py, 0.3.15) — consistent with the rest of the pipeline's `quay.io/biocontainers/*` convention
  - `mgibio/hap.py:v0.3.12` (Docker Hub, updated ~2 years ago)
  Recommend the biocontainers one for consistency with `fastqc`, `multiqc`, `mlst`, etc. already used across your other projects.

---

## 6. COSMIC Cancer Gene Census — registration is on you

Confirmed COSMIC has moved its licensing/registration flow to **cosmickb.org** (redirects from `cancer.sanger.ac.uk/cosmic/license`) — same registration substance (organisational email, free for academic/non-profit use, commercial use requires a paid license), just a domain change worth noting in `docs/data_sources.md` since older tutorials/links will point at the old domain. Current release is v100 (confirmed via `cancer.sanger.ac.uk/cosmic/download/cosmic/v100/mutantcensus`), ahead of the "v99+" the plan cites in §7.

**This step needs you personally:** registration requires an organisational email and manual approval — I cannot create this account. Once registered, download the Cancer Gene Census TSV, note the exact version in `docs/data_sources.md`, and confirm it lands in `.gitignore` (already planned in repo structure, §9 of the plan) before anything else touches it.

---

## 7. Compute / storage estimate (provisional — refine with real ENA numbers)

Without the exact PRJEB27698 Illumina run sizes (blocked by ENA's rate limit this session — see §3 above), here's a standard-WGS-scaling estimate to plan disk headroom against your 12GB-WSL-memory / non-disk-constrained-in-principle setup:

- Typical short-read WGS at ~40–100X combined tumour+normal coverage: **FASTQ ~150–350 GB** (gzipped, both samples combined), **BAM (post-dedup) ~120–280 GB**, **plus ~50–100 GB** for intermediate Mutect2/CNVkit/VEP outputs per full run.
- **Total working-space recommendation for a full-WGS run: 500 GB–700 GB free**, on top of reference resources (GRCh38 FASTA + BWA-MEM2 index ~15 GB, GATK PoN + gnomAD resource VCFs ~15–20 GB combined, VEP cache ~20 GB).
- **Mutect2 runtime at 4 CPU / 8 GB RAM (your WSL2 config):** full-genome tumour-normal Mutect2 typically runs several hours to ~1 day per sample pair at this resource ceiling; this is the concrete reason the plan calls for a `dev` profile before committing to full-WGS turnaround.

**Action for you:** once you have the real per-run base counts from the `curl`/`awk` command in §3, I (or you) can turn this into an exact figure — that's a five-minute follow-up, not a blocker to starting Phase 1 on the `dev` profile.

---

## 8. Nextflow `dev`/`full` profiles — scaffolded

Added to `nextflow.config` (delivered in this batch): a `dev` profile constrained to **chr21 + a melanoma driver gene interval list** (`data/gene_lists/melanoma_genes.tsv` — seeded with the classic COLO829-relevant genes: BRAF, NRAS, CDKN2A, TERT promoter, PTEN, TP53; extend once Module 8's COSMIC CGC integration is live) and a `full` profile with no interval restriction. Both inherit the `executor.cpus = 4` / `executor.memory = '8 GB'` ceiling from your WSL2 config, consistent with your other Nextflow pipelines.

**Not yet done:** the actual interval BED/list file with real GRCh38 coordinates for those gene loci — that's a Phase 1 task once alignment is being built, not a Phase 0 blocker, since profile *structure* was what Phase 0 asked for.

---

## 9. `run_manifest.json` schema — drafted

Delivered as `docs/run_manifest.schema.json` (draft) plus a worked example, covering every field the plan's §9 repo-structure listing named (reference build, Mutect2/VEP/gnomAD/COSMIC CGC/SigProfiler versions, truth set citation + URL, PoN source + URL, pipeline git commit) plus a couple of provenance fields the plan's field list implied but didn't spell out (run timestamp, Nextflow version, container digests) — see that file for the full schema and rationale per field.

---

## Summary — what's actually blocking Phase 1

Nothing here blocks starting Phase 1 (QC/alignment, Modules 1–2) on the `dev` profile against the PRJEB27698 FASTQs. What **should** be resolved before Phase 3 (benchmarking, the project's centrepiece) is the truth-set decision in §1 — that determines what `docs/benchmarking_results.md` is actually measuring against. Recommend deciding that once you've had a chance to look at the options above, rather than defaulting silently.

**Needs your direct action (cannot be done from this sandbox):**
1. Decide the truth-set strategy (§1)
2. Register at COSMIC/cosmickb.org and download the CGC TSV (§6)
3. Run the ENA `curl` command from WSL to get exact FASTQ run accessions and sizes (§3)
4. `docker pull` and confirm the CNVkit and hap.py image tags on your machine (§5)
5. Confirm disk space against the estimate in §7 once real sizes are in hand
