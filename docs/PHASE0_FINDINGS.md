# Phase 0 — Pre-Implementation Findings
**Project 8 — COLO829 Tumour-Normal Somatic Genomics Pipeline**
**Date:** 2026-08-29 (four updates same day: truth-set decision recorded; NYGC's dead download URL tracked down to a working GitHub source and a GRCh37/38 question resolved; the file actually downloaded and extracted successfully; FASTQ run-accession ambiguity found and resolved via independent corroboration)
**Status:** Truth set is decided, downloaded, and extracted. FASTQ run accessions are resolved (`ERR2752449` normal / `ERR2752450` tumour). Nothing left in this document blocks Phase 1 or Phase 3 — only the actual FASTQ download and first `nextflow run` remain, both on your machine.

This document reports what was resolved against the Phase 0 checklist in the v1.1 plan (§3 and §11). Full citations and download commands are in `docs/data_sources.md`. Items requiring your direct action (registration, actual downloads, running commands against your own network) are marked accordingly — none of that can be done from this sandbox.

---

## ✅ 1. Truth set — citation corrected, decision made

The plan names "SEQC2 / Fang et al., *Nature Biotechnology* 2021" as the COLO829 truth set, flagged "to be verified." It doesn't hold up:

**Fang et al. 2021 (the actual SEQC2 Somatic Mutation Working Group flagship paper) benchmarks HCC1395/HCC1395BL, not COLO829.** Its released call sets, at `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/`, are for HCC1395. COLO829 doesn't appear in that paper or on the SEQC2 consortium site (sites.google.com/view/seqc2) as one of the working group's reference pairs. This is exactly the kind of thing Phase 0 exists to catch before it propagates into `docs/benchmarking_results.md`.

The real options for a COLO829 SNV/indel truth set are less clean than the plan assumed:

| Source | Scope | Access | Genome build | Fit |
|---|---|---|---|---|
| **Craig et al. 2016**, *Sci Rep* (["A somatic reference standard for cancer genome sequencing"](https://www.nature.com/articles/srep24607)) | SNV + indel + CNV, 3-platform consensus (>35,000 point mutations, 446 indels) | **Confirmed controlled-access only** — the paper's data availability statement puts "all BAMs and VCFs, including that for the final somatic reference" behind dbGaP (phs000932) or EGA ([EGAS00001001385](https://ega-archive.org/studies/EGAS00001001385), DAC [EGAC00001000408](https://www.omicsdi.org/dataset/ega/EGAC00001000408)). No open supplementary-file shortcut exists. | Not confirmed (2016-era, likely GRCh37 originally) | Widely cited as "the" COLO829 gold standard, but genuinely not usable without a DAC application |
| **Valle-Inclan et al. 2022**, *Cell Genomics* (["A multi-platform reference for somatic structural variation detection"](https://www.sciencedirect.com/science/article/pii/S2666979X22000726)) | **Structural variants only** — no SNV/indel truth set | Fully open, ENA [PRJEB27698](https://www.ebi.ac.uk/ena/browser/view/PRJEB27698), truth VCFs on Zenodo ([4716169](https://zenodo.org/records/4716169), GRCh37 + GRCh38-lifted) | GRCh37 native, GRCh38-lifted VCF also provided | Same sample lineage as the open FASTQs below, but doesn't cover what Module 5 needs (SNV/indel) |
| **NYGC / Xiao et al. 2019**, *Sci Rep* (["Deep whole-genome sequencing of 3 cancer cell lines on 2 sequencing platforms"](https://www.nature.com/articles/s41598-019-55636-3)) | SNV/indel/CNV/SV somatic calls, HiSeqX + NovaSeq, up to 278X | **Processed somatic VCFs open, no login**, via the [NYGC companion site](https://www.nygenome.org/bioinformatics/3-cancer-cell-lines-on-2-sequencers/); raw FASTQ/BAM on **dbGaP phs001839 (controlled access)** | GRCh38 | Open SNV/indel calls, but from NYGC's own pipeline (not an independent multi-platform consensus like Craig 2016) and from different library preps than PRJEB27698 |

**Decision (2026-08-29): NYGC open somatic VCF (HiSeqX version) is the primary SNV/indel truth set.** Reasoning: it's the only genuinely open SNV/indel option, HiSeqX platform-matches the PRJEB27698 FASTQs (removing one axis of difference even though the sequencing centre/library prep still differ), and it keeps the whole project free of controlled-access dependencies — preserving the exact rationale that motivated choosing COLO829 over TCGA in the first place. Pursuing EGA DAC access for Craig et al. 2016 was explicitly declined for that reason.

The lineage-mismatch caveat (NYGC's own pipeline, separate sequencing run, not an independent multi-platform consensus like Craig 2016) is now recorded in `docs/data_sources.md` §3, to carry into `docs/benchmarking_results.md` per the plan's Phase 0 instruction not to leave this quietly assumed. Valle-Inclan et al. 2022's SV truth set (fully open, genuinely PRJEB27698-lineage-matched) is kept as a complementary Module 7 CNV/SV cross-check, not a substitute — it doesn't cover SNV/indel so it can't replace the NYGC VCF for Module 5.

**Update 2026-08-29 (later same day): the NYGC page in the citation above is dead.** The paper cites `www.nygenome.org/bioinformatics/...`; NYGC has since moved their bioinformatics site to a `bioinformatics.nygenome.org` subdomain with no redirect from the old URL, and even the new subdomain's own download buttons didn't resolve through automated fetching (possibly JS-rendered rather than plain links) or, per your own browser test, sometimes fail outright. Tracked down a working alternative instead: **NYGC's own GitHub repo**:
```
https://github.com/nygenome/3-cancer-cell-lines-on-2-sequencers
https://raw.githubusercontent.com/nygenome/3-cancer-cell-lines-on-2-sequencers/master/data/Variants.HighCoverage.tar.gz
```

**Update 2026-08-29 (later still): downloaded and extracted successfully.** One snag along the way, resolved: the archive downloaded correctly (73,643,008 bytes, matching GitHub's `content-length` header exactly — not truncated, not an HTML error page), but `tar xzf` failed with "not in gzip format." Turned out the file is a **plain (non-gzip) tar archive despite the `.tar.gz` name** — README's "140MB" was also stale, real size is ~70.2MB. Fix was simply `tar xf` instead of `tar xzf` (GNU tar auto-detects format). Confirmed file now in hand: `COLO-829--COLO-829BL.snv.indel.final.v6.annotated.vcf` (36,165,120 bytes, uncompressed). Also confirmed while extracting: NYGC does **not** publish a high-confidence callable-regions BED for SNV/indel (only CNV bed + SV bedpe exist per cell line) — hap.py will run without a confidence-region restriction unless a callability filter is added later, a Phase 3 decision. Exact commands and file paths are in `docs/data_sources.md` §3.

**A GRCh37-vs-GRCh38 question came up while chasing this down and was worth checking properly rather than taking at face value:** the paper separately re-ran its pipeline on GRCh37, purely to compare against Craig et al. 2016's GRCh37-only, EGA-controlled dataset (EGAD00001002142 — the same Craig 2016 study already declined above, not a new resource). Confirmed by reading the paper's actual methods text: **the primary pipeline output we're downloading — the GRCh38 `Variants.HighCoverage.tar.gz` file — is unaffected**; the GRCh37 run is a paragraph in their methods, not a file we touch. If anything this strengthens the choice: the paper reports 98% concordance between NYGC's own calls and the Craig et al. 2016 gold-standard set, i.e. an independent validation of the exact file we're using, already done for us. Full detail in `docs/data_sources.md` §3.

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

## 3. FASTQ accessions — resolved (2026-08-29, second pass)

**ENA project [PRJEB27698](https://www.ebi.ac.uk/ena/browser/view/PRJEB27698)** ([Valle-Inclan et al. 2022](https://www.biorxiv.org/content/10.1101/2020.10.15.340497v1.full) data availability statement) is confirmed **fully open access, no approval required**, and includes Illumina HiSeq X Ten short-read WGS for both COLO829 (tumour) and COLO829BL (normal) — this is the same project that produced the open SV truth set above, so FASTQ and (SV) truth set share lineage, which is a plus.

**Caveat:** PRJEB27698 is a *multi-platform* project — 29 SRA experiments, ~4.09 Tbases / 1.55 TB total, spanning Illumina HiSeq, 10x Genomics linked-reads, Oxford Nanopore, PacBio, and BioNano optical mapping. **Do not bulk-download the whole project** — only the plain Illumina HiSeq X Ten WGS runs are relevant to this pipeline's BWA-MEM2/Mutect2 design.

**First pass (blocked from this sandbox):** ENA's API rate-limited this sandbox's WebFetch repeatedly, so the run-accession-level query had to be run from your machine instead:
```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB27698&result=read_run&fields=run_accession,sample_title,instrument_platform,library_strategy,fastq_ftp,base_count&format=tsv" \
  | awk -F'\t' '$3=="ILLUMINA"'
```

**Second pass — ambiguity found and resolved.** Your query returned 9 ILLUMINA/WGS rows. Five (`ERR4093255`–`ERR4093259`) are a purity-titration series, out of scope. The remaining four split into two same-label pairs — "COLO829 reference cell line" (`ERR2752449` vs. `ERR2820166`) and "COLO829 melanoma cell line" (`ERR2752450` vs. `ERR2820167`) — with no obvious way to tell which of each pair is the real plain WGS run just from the ENA row itself. Rather than guess (the same mistake risk that bit the truth-set citation in §1), this was resolved via an independent source: [Cameron et al. 2021, "GRIDSS2"](https://www.biorxiv.org/content/10.1101/2020.07.09.196527v2.full) (*Genome Biology*, same dataset/collaboration as Valle-Inclan et al. 2022) states in its methods that it used *"Illumina HiSeqX (ENA run accessions ERR2752449 and ERR2752450 for COLO829BL and COLO829T, respectively)."* This also lines up with the fact that only that pair has a populated `fastq_ftp` field — the other pair (`ERR2820166`/`ERR2820167`) is most likely a suppressed/superseded duplicate submission.

**Decision: `ERR2752449` = COLO829BL (normal), `ERR2752450` = COLO829 (tumour).** Full accession table, download commands, and implied per-run coverage (~37X normal / ~98X tumour): `docs/data_sources.md` §1.

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

**Update 2026-08-29: registration attempted, deferred — not currently blocking.** COSMIC's account form asks for an organisation type (For-Profit Company / Not-For-Profit Hospital / Not-For-Profit Institution-Academic Research Centre / Private Hospital); this project is a personal, non-commercial portfolio effort with no institutional affiliation, so none of those categories cleanly fits a private individual, and registration is on hold for now. This is fine to leave open: **COSMIC's Cancer Gene Census is only consumed in Phase 5** (the oncogenicity/actionability interpretation layer, cross-referencing CGC driver genes against called variants) — Phases 1-4 (QC/alignment, contamination + Mutect2, benchmarking, CNVkit) never touch it. If registration is still unresolved when Phase 5 starts, [CIViC](https://civicdb.org) is a fully open, no-registration alternative covering similar oncogenicity/actionability ground — not a like-for-like replacement for CGC's driver-gene census, but a legitimate documented fallback rather than a hard stop. If a paying client engagement materialises before then, registering as "For-Profit Company" with COSMIC's paid commercial license at that point is the honest path, not something to route around.

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

Nothing blocks starting Phase 1 (QC/alignment, Modules 1–2) on the `dev` profile against the PRJEB27698 FASTQs. The truth set (§1) is fully in hand now — decided, downloaded, and extracted — so Phase 3 (benchmarking) is unblocked too.

**Needs your direct action (cannot be done from this sandbox):**
1. ~~Decide the truth-set strategy~~ — done (§1): NYGC open VCF, HiSeqX version
2. ~~Find a working download URL and download it~~ — done (§1): GitHub source, downloaded and extracted; `COLO-829--COLO-829BL.snv.indel.final.v6.annotated.vcf` confirmed in hand
3. Register at COSMIC/cosmickb.org and download the CGC TSV (§6) — **deferred, not blocking**: not needed until Phase 5 (oncogenicity/actionability interpretation); registration is on hold pending either a clearer fit for a private-individual account or a future commercial engagement, with CIViC noted as an open fallback if it's still unresolved by then
4. ~~Run the ENA `curl` command from WSL to get exact FASTQ run accessions and sizes~~ — done (§3): `ERR2752449` (normal, ~114 Gbases) / `ERR2752450` (tumour, ~303 Gbases); download commands in `docs/data_sources.md` §1
5. `docker pull` and confirm the CNVkit and hap.py image tags on your machine (§5) — `samtools:1.21--h50ea8bc_0` already confirmed working
6. Confirm disk space against the ~210–240 GB combined raw-FASTQ + ~3 GB reference estimate in `docs/data_sources.md` §1/§2 (refines the older §7 rule-of-thumb now that real base counts are in hand)
7. ~~Resolve which GRCh38 reference FASTA to use~~ — done (`docs/data_sources.md` §2): `Homo_sapiens_assembly38.fasta`, confirmed to be the correct pairing for the PoN/gnomAD files already in `nextflow.config`. Download links given there are **unverified from this sandbox** — run the `--list-only` FTP check first before pulling the full ~3 GB file.
8. Actually download the two FASTQ pairs and the reference (`docs/data_sources.md` §1/§2 download commands), generate `.fai`/`.dict`, and run the first real `nextflow run main.nf -profile docker,dev ...` per `docs/PHASE1_NOTES.md`
