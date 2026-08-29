# Data Sources, Versions, and Licensing

This file is the single source of truth for every external resource this pipeline depends on: where it comes from, what version was actually used, and what its access/licensing terms are. Populate the "confirmed" columns as each resource is actually downloaded (this is a living document, not a one-time Phase 0 artifact) — the goal is that `results/provenance/run_manifest.json` (see `docs/run_manifest.schema.json`) and this file always agree.

Full research trail and rationale for the decisions below: `docs/PHASE0_FINDINGS.md`.

---

## 1. Sequencing data — COLO829 / COLO829BL FASTQs

- **Source:** ENA project [PRJEB27698](https://www.ebi.ac.uk/ena/browser/view/PRJEB27698)
- **Citation:** Valle-Inclan, J.E. et al. "A multi-platform reference for somatic structural variation detection." *Cell Genomics* 2022. https://www.sciencedirect.com/science/article/pii/S2666979X22000726 (preprint: https://www.biorxiv.org/content/10.1101/2020.10.15.340497v1.full)
- **Access:** Open, no approval required (confirmed in the paper's data availability statement)
- **Platform needed:** Illumina HiSeq X Ten short-read WGS only — PRJEB27698 is a multi-platform project (29 SRA experiments, ~4.09 Tbases total across Illumina, 10x Genomics, Oxford Nanopore, PacBio, BioNano). **Do not bulk-download the whole project.**
- **Run accessions used:** TODO — fill in after running the ENA filereport query below from a machine with unrestricted network access (this sandbox's WebFetch was rate-limited by ENA repeatedly):
  ```bash
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB27698&result=read_run&fields=run_accession,sample_title,instrument_platform,library_strategy,fastq_ftp,base_count&format=tsv" \
    | awk -F'\t' 'NR==1 || $3=="ILLUMINA"'
  ```
- **Coverage / instrument / library prep:** TODO — record here once the query above is run
- **File sizes (download + expected BAM):** TODO — see `docs/PHASE0_FINDINGS.md` §7 for a provisional rule-of-thumb estimate pending real numbers

---

## 2. Reference genome

- **Build:** GRCh38 (specific decoy/patch level TODO — confirm which GRCh38 variant, e.g. `GRCh38.d1.vd1` vs. plain primary assembly, matches what the GATK PoN/gnomAD resources below were built against, to avoid contig-mismatch errors in Mutect2)
- **Source:** TODO — standard Broad/GATK GRCh38 reference bundle, download command + md5 to be recorded here in Phase 1

---

## 3. Somatic truth set — **DECIDED and DOWNLOADED 2026-08-29**

The plan's original citation (SEQC2 / Fang et al. 2021, *Nature Biotechnology*) is **incorrect for COLO829** — that paper's reference pair is HCC1395/HCC1395BL.

**Decision: NYGC open somatic VCF (HiSeqX version) is the primary SNV/indel truth set for Module 5 benchmarking.**

- **Citation:** Xiao, W. et al. "Deep whole-genome sequencing of 3 cancer cell lines on 2 sequencing platforms." *Sci Rep* 2019. https://www.nature.com/articles/s41598-019-55636-3
- **Data / working download (verified 2026-08-29):** the paper's cited URL (`www.nygenome.org/bioinformatics/...`) is dead — NYGC moved their bioinformatics site to a `bioinformatics.nygenome.org` subdomain and the old link was never redirected. The reliable current source is NYGC's own GitHub repo, confirmed live:
  ```
  https://github.com/nygenome/3-cancer-cell-lines-on-2-sequencers
  ```
  Direct download and extraction, **actually run and confirmed working 2026-08-29**:
  ```bash
  mkdir -p ~/projects/somatic-variant-analysis-COLO829/truth_set && cd $_
  curl -LO https://raw.githubusercontent.com/nygenome/3-cancer-cell-lines-on-2-sequencers/master/data/Variants.HighCoverage.tar.gz
  tar tf Variants.HighCoverage.tar.gz | grep -i colo   # sanity check before extracting everything
  tar xf Variants.HighCoverage.tar.gz                  # NOTE: despite the .tar.gz name this is a PLAIN (non-gzip) tar — do not use -z, GNU tar auto-detects
  ```
  Extract only the **COLO-829 (HiSeqX)** SNV/indel VCF + CNV BED, not the NovaSeq version: HiSeqX is a platform match to the PRJEB27698 FASTQs (also HiSeq X Ten), which removes one axis of difference even though sequencing centre and library prep still differ.
- **Confirmed file in hand:** `~/projects/somatic-variant-analysis-COLO829/truth_set/COLO-829--COLO-829BL.snv.indel.final.v6.annotated.vcf` — 36,165,120 bytes, uncompressed VCF (not `.vcf.gz`). For hap.py you'll likely want it bgzip-compressed and tabix-indexed:
  ```bash
  bgzip -k COLO-829--COLO-829BL.snv.indel.final.v6.annotated.vcf
  tabix -p vcf COLO-829--COLO-829BL.snv.indel.final.v6.annotated.vcf.gz
  ```
- **Companion CNV file also extracted:** `COLO-829--COLO-829BL.cnv.annotated.v6.final.bed` (13,490 bytes) — not the truth set for Module 5, but usable later as an extra Module 7 cross-check alongside Valle-Inclan's SV truth set.
- **Access:** Open, no login (GitHub-hosted, no dbGaP/EGA involved for this file)
- **Genome build:** **GRCh38**, confirmed from the paper's methods text — "Sequencing reads were aligned to the GRCh38 reference genome (1000 Genomes version)." (Exact decoy/patch level still TODO — see §2 below, must match what Mutect2/PoN/gnomAD are built against.)
- **Important clarification (resolved 2026-08-29):** the paper *separately* re-ran their pipeline on GRCh37 purely to compare against Craig et al. 2016's GRCh37-only, EGA-controlled-access truth set (EGAD00001002142, part of EGAS00001001385 — the same Craig 2016 dataset already declined below, not a new resource). That GRCh37 comparison run is not a file we use; it's a paragraph in their methods. **The GRCh38 `Variants.HighCoverage.tar.gz` file above is what we download**, and it's unaffected by the GRCh37 side-comparison. That comparison is actually a positive signal for our choice: the paper reports "we called over 98% of the Craig et al. SNVs" — i.e. NYGC's own calls were independently validated at 98% concordance against the real gold-standard Craig 2016 set, which is exactly the multi-platform consensus truth set we can't access directly (controlled-access).
- **Bundle:** `Variants.HighCoverage.tar.gz` — actual size **73,643,008 bytes (~70.2 MiB)**, confirmed via HTTP `content-length` header matching the downloaded file exactly (README's stated "140MB" was stale/wrong, not a download problem). Contains all three cell lines' HiSeqX + NovaSeq VCFs, CNV BEDs, and SV bedpe files — only the COLO-829 HiSeqX files are relevant here.
- **High-confidence BED — resolved:** NYGC does **not** publish a callable-regions BED for SNV/indel — the archive only contains a CNV bed and SV/SV-high-confidence bedpe per cell line, no SNV/indel-specific confidence regions. hap.py benchmarking will run without a confidence-region restriction unless a separate callability/mappability filter is applied — a Phase 3 decision, not a Phase 0 blocker.
- **Raw NYGC sequencing data (for reference, not used by us):** dbGaP study `phs001839.v1.p1` — https://dbgap.ncbi.nlm.nih.gov/beta/study/phs001839.v1.p1/ — controlled access, irrelevant to this pipeline since we source FASTQs from PRJEB27698 instead (§1 above)

**Caveat, carried into `docs/benchmarking_results.md`:** this truth set was generated by NYGC's own pipeline from a separate HiSeqX sequencing run/library prep, not the exact PRJEB27698 FASTQs this pipeline aligns. It is not itself an independent multi-platform consensus (Craig et al. 2016 is that, but is EGA/dbGaP controlled-access only) — it's one lab's own somatic calls. However, per the clarification above, NYGC's own paper already cross-validated these calls at 98% concordance against Craig et al. 2016, which meaningfully strengthens confidence in using it as our benchmark. Still, frame results against it as "agreement with a GATK-family-adjacent somatic pipeline, itself independently validated against the field's gold standard" rather than as unmediated ground truth — the exact caveat the plan's Phase 0 checklist asked to be documented, not quietly assumed.

**Why not Craig et al. 2016** (Craig, D.W. et al., *Sci Rep* 2016, https://www.nature.com/articles/srep24607): confirmed 2026-08-29 that the paper's own data availability statement puts *all* BAMs and VCFs — "including that for the final somatic reference" — behind dbGaP (phs000932) or EGA (EGAS00001001385) controlled access, with no open supplementary-file shortcut. Despite being the most-cited 3-platform consensus COLO829 truth set, using it would require a formal Data Access Committee application and would break the project's "no controlled-access issues" framing that motivated choosing COLO829 over TCGA in the first place. Deliberately not pursued.

**Complementary, not a substitute — Valle-Inclan et al. 2022 SV truth set:** fully open (see §1 above; truth VCF at https://zenodo.org/records/4716169), same PRJEB27698 lineage as the input FASTQs. Doesn't cover SNV/indel, so it doesn't replace the NYGC VCF for Module 5 — but it's a good open cross-check for Module 7's CNVkit output, since it's genuinely lineage-matched in a way the NYGC calls aren't. Worth wiring in as an optional Module 7 comparison once CNVkit is running.

---

## 4. GATK somatic resources (hg38)

Source: Broad public GCS bucket `gs://gatk-best-practices/somatic-hg38/` (also mirrored at `https://storage.googleapis.com/gatk-best-practices/somatic-hg38/<filename>` for non-GCP-authenticated download).

| Resource | Path | Used for |
|---|---|---|
| Panel of Normals | `gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz` | Mutect2 `--panel-of-normals` |
| gnomAD germline resource | `gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz` | Mutect2 `--germline-resource` |
| Common biallelic sites | `gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz` | GetPileupSummaries |

**gnomAD version actually embedded in `af-only-gnomad.hg38.vcf.gz`:** TODO — check the VCF header after download; the plan's §8 cites "gnomAD v4.1" but the Broad bucket file's actual embedded version needs confirming, they are not necessarily the same curation.

---

## 5. COSMIC Cancer Gene Census

- **Registration:** https://cancer.sanger.ac.uk/cosmic/register (redirects to cosmickb.org's licensing/registration flow — domain changed since the plan was written)
- **License:** Free for academic/non-profit use with an organisational-email account; commercial use requires a separate paid license. Full terms: https://www.cosmickb.org/terms/
- **Current release at time of writing:** v100 (plan's §7 cites "v99+" — update once you've actually registered and downloaded, record the exact version pulled)
- **Repository handling:** CGC TSV is **never committed** — `.gitignore` excludes `COSMIC_*.tsv`; only parsing code and these download instructions live in the repo (mirrors the PGx project's CPIC/PharmVar handling)
- **Downloaded version / date:** TODO — fill in once you've registered and pulled the file

---

## 6. Mutational signature tools — corrected package names

The plan's Module 9 spec names `SigProfilerSingleBase`, which **does not exist on PyPI** (confirmed by a failed `pip install` in this session). The correct Alexandrov Lab packages for "decompose the observed 96-channel spectrum against COSMIC SBS reference set v3.4":

- `SigProfilerMatrixGenerator` (PyPI, latest 1.3.6) — builds the 96-channel SBS matrix from the PASS VCF
- `SigProfilerAssignment` (PyPI, latest 1.1.5) — attributes the observed spectrum against a reference signature set (this replaces the "attribution" role the plan assigned to the nonexistent package)

Pin exact versions in `requirements.txt` once Module 9 is implemented (Phase 6).

---

## 7. CNVkit / hap.py containers

| Tool | Image | Status |
|---|---|---|
| CNVkit | `quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0` | Plan's original choice confirmed valid (0.9.10-0 exists in bioconda); verify exact tag with a live `docker pull` before wiring into `cnvkit.nf`, since biocontainers build-hash suffixes can shift between rebuilds |
| hap.py | `quay.io/biocontainers/hap.py:0.3.15-0` (recommended, replaces the plan's `pkrusche/hap.py:latest`) | `pkrusche/hap.py:latest` confirmed unmaintained (~9 years since last push) — switched to the actively-rebuilt biocontainers image for consistency with the rest of the pipeline's `quay.io/biocontainers/*` convention |

---

## 8. BQSR decision (carried over from the plan, unchanged)

Base Quality Score Recalibration is **deliberately not applied** before Mutect2. Per the Broad's current somatic best-practices guidance, BQSR's benefit for somatic calling is less clear than for germline, and is considered optional for WGS cell-line data without FFPE artefacts. Documented here per the plan's instruction not to leave this quietly assumed.

---

## 9. gnomAD population AF resource for VEP annotation

Same open question as item 4 (gnomAD v4.1 vs. whatever's embedded in the Broad bucket file) — Module 6's annotation step and Module 4's germline-resource filtering should cite a consistent gnomAD version; confirm and record here once both are downloaded.
