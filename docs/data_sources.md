# Data Sources, Versions, and Licensing

This file is the single source of truth for every external resource this pipeline depends on: where it comes from, what version was actually used, and what its access/licensing terms are. Populate the "confirmed" columns as each resource is actually downloaded (this is a living document, not a one-time Phase 0 artifact) — the goal is that `results/provenance/run_manifest.json` (see `docs/run_manifest.schema.json`) and this file always agree.

Full research trail and rationale for the decisions below: `docs/PHASE0_FINDINGS.md`.

---

## 1. Sequencing data — COLO829 / COLO829BL FASTQs

- **Source:** ENA project [PRJEB27698](https://www.ebi.ac.uk/ena/browser/view/PRJEB27698)
- **Citation:** Valle-Inclan, J.E. et al. "A multi-platform reference for somatic structural variation detection." *Cell Genomics* 2022. https://www.sciencedirect.com/science/article/pii/S2666979X22000726 (preprint: https://www.biorxiv.org/content/10.1101/2020.10.15.340497v1.full)
- **Access:** Open, no approval required (confirmed in the paper's data availability statement)
- **Platform needed:** Illumina HiSeq X Ten short-read WGS only — PRJEB27698 is a multi-platform project (29 SRA experiments, ~4.09 Tbases total across Illumina, 10x Genomics, Oxford Nanopore, PacBio, BioNano). **Do not bulk-download the whole project.**
- **Run accessions used — RESOLVED 2026-08-29:**
  Your ENA `curl`/`awk` query returned 9 ILLUMINA/WGS rows under PRJEB27698. Five of those (`ERR4093255`–`ERR4093259`) are a purity-titration dilution series (10/20/25/50/75% tumour purity, ~300–304 Gbases each) — **out of scope**, set aside. The remaining four split into two ambiguous same-label pairs:

  | Label | Accession | Base count | `fastq_ftp` |
  |---|---|---|---|
  | "COLO829 reference cell line" (normal) | `ERR2752449` | 114,093,383,622 | present |
  | "COLO829 reference cell line" (normal) | `ERR2820166` | 99,480,531,208 | **empty** |
  | "COLO829 melanoma cell line" (tumour) | `ERR2752450` | 303,329,725,932 | present |
  | "COLO829 melanoma cell line" (tumour) | `ERR2820167` | 103,505,548,977 | **empty** |

  This was deliberately **not guessed** — instead resolved by independent corroboration: [Cameron et al. 2021, "GRIDSS2"](https://www.biorxiv.org/content/10.1101/2020.07.09.196527v2.full) (*Genome Biology*, same PRJEB27698 dataset, same research collaboration as Valle-Inclan et al. 2022) states explicitly in its methods: *"Illumina HiSeqX (ENA run accessions ERR2752449 and ERR2752450 for COLO829BL and COLO829T, respectively)."*

  **Decision: `ERR2752449` = COLO829BL (normal), `ERR2752450` = COLO829 (tumour).** This also matches the internally-consistent signal that only this pair has a populated `fastq_ftp` field — the `ERR2820166`/`ERR2820167` pair is most likely a superseded or otherwise-suppressed duplicate submission, not a second legitimate dataset. Implied coverage (base_count ÷ 3.1 Gb haploid genome): normal ≈ **36.8X**, tumour ≈ **97.8X** — a plausible, somewhat-higher-tumour-depth WGS design for somatic sensitivity, consistent with the paper's WGS design intent (not the same number as the paper's own quoted "155X"/"235X" combined-coverage figures, which sum across all five sequencing technologies in the project, not just this one HiSeqX run — see caveat below).
- **Download commands** (run from WSL — get the exact `fastq_ftp` URLs for just these two accessions, don't bulk-download the project):
  ```bash
  mkdir -p ~/projects/somatic-variant-analysis-COLO829/fastq && cd $_
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB27698&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,base_count&format=tsv" \
    | awk -F'\t' '$1=="ERR2752449" || $1=="ERR2752450" || NR==1'
  # fastq_ftp column is two semicolon-separated URLs (R1;R2) per run -- prefix each with https:// and curl -LO them,
  # or feed straight to wget -c. Verify against the fastq_md5 column after download.
  ```
- **Coverage / instrument / library prep:** Illumina HiSeq X Ten, TruSeq Nano prep (per Valle-Inclan et al. 2022 methods, §3 below). Per-run depth: see table above.
- **File sizes (download + expected BAM):** normal ~114 Gbases / tumour ~303 Gbases of raw read data — refines the provisional rule-of-thumb estimate in `docs/PHASE0_FINDINGS.md` §7 into something concrete: expect roughly 60–70 GB gzipped FASTQ for the normal pair and 150–170 GB for the tumour pair (typical ~0.5–0.6 bytes/base for gzipped Illumina short reads), i.e. **~210–240 GB combined raw FASTQ**, before BAM/dedup/Mutect2 intermediates. Confirm against actual downloaded file sizes once fetched.
- **Caveat carried into `docs/benchmarking_results.md`:** the paper's own reported "235X"/"155X" (tumour/normal) coverage figures are *combined across all five sequencing technologies* used in PRJEB27698 (Illumina, 10x Genomics, Nanopore, PacBio, BioNano) — not this HiSeqX run alone. Don't quote those numbers as this pipeline's actual input depth; use the per-run `base_count`-derived ~37X/~98X figures above, and replace with `samtools depth`/`stats` output once BAMs exist (Phase 1).
- **Donor sex — confirmed 2026-08-31 (needed for Phase 4/CNVkit's sex-chromosome ploidy handling):** [ATCC's own COLO 829 (CRL-1974) product page](https://www.atcc.org/products/crl-1974) lists the donor as a 45-year-old **male**. Verified rather than assumed, since getting this wrong would silently bias every CNVkit chrX/chrY call — `modules/cnvkit.nf` passes `-y`/`--male-reference` on that basis. See `docs/PHASE4_NOTES.md`.

**Real-depth data via ENA's pre-aligned GRCh37 BAMs — found and extracted 2026-09-01, see `docs/MUTECT2_SCATTERGATHER_NOTES.md`'s context for why this was needed:** every dedup BAM this project had actually produced through 2026-09-01 turned out to be the same 10,000-read-pair `dev`-profile subsample, confirmed via `samtools flagstat` — despite `docs/PHASE1_NOTES.md`/the README describing Phase 1 as having "run successfully end-to-end against the real full genome." That phrase meant the full GRCh38 *reference* was used as the alignment target, not that real full-sequencing-depth reads existed — a real ambiguity in that earlier wording that led to a wrong assumption later in the project (see `docs/MUTECT2_SCATTERGATHER_NOTES.md`). No real full-depth FASTQ or BAM had ever actually been downloaded.

Checking ENA's file listing for the real accessions (§1 above) found something better than the raw FASTQ pull originally planned: both runs have a **pre-aligned, deduplicated, realigned BAM submitted directly to ENA**, each with a `.bai` index —
```
ERR2752450 (tumour, COLO829T): https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752450/COLO829T_dedup.realigned.bam(.bai)
ERR2752449 (normal, COLO829R): https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752449/COLO829R_dedup.realigned.bam(.bai)
```
confirmed via `curl` against ENA's own filereport API (`submitted_ftp`/`submitted_format` columns), not assumed. With a `.bai` present, `samtools` can pull just the reads overlapping a chosen region straight over HTTPS (range requests against the remote index) — no need to download or align the full ~174GB combined raw FASTQ (exact sizes confirmed via the same API call's `fastq_bytes` column: ~124.6GB tumour, ~49.4GB normal) for a first real-depth test.

**Reference build mismatch — checked, not assumed:** `samtools view -H` against the tumour BAM shows bare contig names (`SN:1`, `SN:2`, ... `SN:X`/`SN:Y`/`SN:MT`) with `LN:249250621` for contig `1` — GRCh37's chr1 length (248,956,422 for GRCh38) confirms these BAMs are aligned to **GRCh37**, not the GRCh38 build (§2 above) the rest of this pipeline uses. Using them directly against `Homo_sapiens_assembly38.fasta`/the GATK hg38 resource bundle would hit exactly the "incompatible contigs" problem §2's own reference-choice caveat warns about, just from the opposite direction (wrong build, not wrong hg38 variant).

**Resolution: use the remote BAM only to extract which real reads land in the melanoma gene panel, then re-align those reads to GRCh38 through this pipeline's own Module 1-2 — not use the BAM's alignment or coordinates directly.** `data/gene_lists/dev_intervals_grch37.bed` (new file, same 8 genes as `dev_intervals.bed`, same Ensembl-lookup-then-pad method but against Ensembl's GRCh37 archive API and confirmed `assembly_name: "GRCh37"` per gene) gives `samtools view -L` the right coordinates for these specific remote BAMs.

**Original extraction commands (single-shot, whole-BED `-L`) — SUPERSEDED 2026-09-02, do not use, kept below only for the paper trail:**
```bash
mkdir -p ~/projects/somatic-variant-analysis-COLO829/real_data && cd $_

# Tumour (ERR2752450 / COLO829T) -- real ~98X depth in-panel
samtools view -b -L ../data/gene_lists/dev_intervals_grch37.bed \
  https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752450/COLO829T_dedup.realigned.bam \
  -o COLO829_real_genepanel.bam
samtools sort -n -o COLO829_real_genepanel.namesorted.bam COLO829_real_genepanel.bam
samtools fastq -1 COLO829_real_R1.fastq.gz -2 COLO829_real_R2.fastq.gz -0 /dev/null -s /dev/null -n COLO829_real_genepanel.namesorted.bam

# Normal (ERR2752449 / COLO829R) -- real ~37X depth in-panel
samtools view -b -L ../data/gene_lists/dev_intervals_grch37.bed \
  https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752449/COLO829R_dedup.realigned.bam \
  -o COLO829BL_real_genepanel.bam
samtools sort -n -o COLO829BL_real_genepanel.namesorted.bam COLO829BL_real_genepanel.bam
samtools fastq -1 COLO829BL_real_R1.fastq.gz -2 COLO829BL_real_R2.fastq.gz -0 /dev/null -s /dev/null -n COLO829BL_real_genepanel.namesorted.bam
```

**Why superseded — overnight run silently truncated (found 2026-09-02):** the single-shot commands above were left running unattended overnight against both remote BAMs. They "completed" (exit 0, files written) but stderr showed `[E::bgzf_read_block] Failed to read BGZF block data` / `Broken pipe` partway through, and — more importantly, since a completed command isn't the same as a correct result — the resulting read counts were not just low but disproportionate between samples: tumour came back with 456,807 reads, normal with only 4,891. Real in-panel depth is ~98X tumour vs ~37X normal, so a clean extraction should land close to a ~2.6:1 ratio; the actual ratio was ~93:1, meaning the normal-sample fetch died almost immediately while the tumour fetch ran much further before also dying. The most likely cause is a dropped connection sometime during an unattended multi-hour single HTTPS stream (WSL/laptop sleep, wifi hiccup, or a server-side idle timeout on `ftp.sra.ebi.ac.uk`) — not a flaw in the region-extraction approach itself, which is otherwise sound (per the boundary-pair caveat below).

**Corrected extraction method: `bin/extract_real_gene_panel.sh`** — splits the single whole-BED fetch into one `samtools view` call per gene (8 small, independently-retried fetches per sample instead of 1 long one), each validated with `samtools quickcheck` before being trusted, with up to 5 retries and a 15s backoff per gene. A single dropped connection now costs one gene's retry, not the whole overnight run, and a truncated/incomplete BGZF stream is caught immediately rather than silently merged into the final FASTQ the way the first attempt's partial output would have been. Usage:
```bash
cd ~/projects/somatic-variant-analysis-COLO829
bash bin/extract_real_gene_panel.sh COLO829   https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752450/COLO829T_dedup.realigned.bam data/gene_lists/dev_intervals_grch37.bed real_data
bash bin/extract_real_gene_panel.sh COLO829BL https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752449/COLO829R_dedup.realigned.bam data/gene_lists/dev_intervals_grch37.bed real_data
```
This writes the same final filenames the original commands targeted (`real_data/COLO829_real_R1/R2.fastq.gz`, `real_data/COLO829BL_real_R1/R2.fastq.gz`), so it overwrites the first attempt's incomplete output directly — no manual cleanup required, though the stray `real_data/*_genepanel*.bam` intermediates from the failed run can be deleted for tidiness. After each sample completes, sanity-check the read-count ratio before trusting the result and moving on to alignment — the whole reason this fix exists is that "it finished" was not sufficient evidence the first time.

The resulting FASTQs feed this pipeline exactly like any other input (`--tumour_reads_1`/etc.), on the `dev` profile (its `-L` restriction to `dev_intervals.bed` — the GRCh38 version — already matches where this data actually has coverage, since reads outside the panel were never extracted).

**Known, disclosed limitation, not silently glossed over:** `samtools view -L`/per-region `samtools view` only requires a read's own alignment position to overlap the region — a read pair straddling a window boundary can have its mate fall just outside it, so a small number of boundary-proximal pairs may come through incomplete (and `samtools fastq`'s `-0`/`-s` disposal of non-properly-paired reads means those get dropped rather than miscounted). The existing ±2000bp padding provides some slack; this isn't treated as a correctness bug, just an honestly-stated edge effect of targeted extraction versus a true full-BAM pull.

---

## 2. Reference genome — RESOLVED 2026-08-29

- **File:** `Homo_sapiens_assembly38.fasta` — the Broad's own GATK hg38 resource-bundle reference (**not** `GRCh38.d1.vd1`, the GDC/TCGA variant with a different decoy set — those two are not interchangeable and mixing them with the wrong PoN/gnomAD file causes exactly the "incompatible contigs" Mutect2 errors this section exists to avoid).
- **Why this file specifically (not just "some GRCh38"):** confirmed via [nf-core/sarek's `igenomes.config`](https://github.com/nf-core/sarek/blob/master/conf/igenomes.config) — an actively-maintained, widely-used GATK-based somatic pipeline — that `Homo_sapiens_assembly38.fasta` is the reference paired with `af-only-gnomad.hg38.vcf.gz` and `1000g_pon.hg38.vcf.gz` under the same GATK bundle "Annotation/GATKBundle" grouping, i.e. the same three filenames this pipeline already uses for `germline_resource` and `panel_of_normals` (§4 below). Same-bundle filenames pairing up across an independent pipeline's config is good evidence these are the versions meant to be used together.
- **Confirmed contents:** includes ALT and decoy contigs (per [Sarek's REFERENCES.md](https://github.com/SciLifeLab/Sarek/blob/master/docs/REFERENCES.md), which lists the accompanying `.alt`/`.amb`/`.ann`/`.bwt`/`.pac`/`.sa` index files) — chr-prefixed naming as used throughout the rest of this pipeline's resources.
- **MD5 (per Sarek's REFERENCES.md, cross-check after download):** `7ff134953dcca8c8997453bbb80b6b5e`
- **Download sources — verify reachability before pulling (~3 GB), this sandbox could not reach any of these to test them directly (network policy blocks `storage.googleapis.com` and FTP from here):**
  ```bash
  # Anonymous FTP (documented in multiple GATK tutorials/Sarek docs) -- try a directory listing first:
  curl -u gsapubftp-anonymous: --list-only ftp://ftp.broadinstitute.org/bundle/hg38/

  # If that lists Homo_sapiens_assembly38.fasta, pull it + companion .fai/.dict:
  curl -u gsapubftp-anonymous: -O ftp://ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta
  curl -u gsapubftp-anonymous: -O ftp://ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai
  curl -u gsapubftp-anonymous: -O ftp://ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict

  # HTTPS mirrors (untested from this sandbox -- try if FTP is blocked on your network instead):
  # https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
  # https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta

  md5sum Homo_sapiens_assembly38.fasta   # compare against 7ff134953dcca8c8997453bbb80b6b5e above
  ```
- **`.fai`/`.dict`:** pull alongside if the FTP directory has them (commands above); if not, generate locally — `samtools faidx Homo_sapiens_assembly38.fasta` and `gatk CreateSequenceDictionary -R Homo_sapiens_assembly38.fasta` — needed by Phase 2's GATK steps regardless (see `docs/PHASE1_NOTES.md` item 5).

---

## 3. Somatic truth set — **DECIDED 2026-08-29, REVISED 2026-08-30 (HiSeqX file found corrupted at source)**

The plan's original citation (SEQC2 / Fang et al. 2021, *Nature Biotechnology*) is **incorrect for COLO829** — that paper's reference pair is HCC1395/HCC1395BL.

**UPDATE 2026-08-30 — the HiSeqX file this section originally chose is unusable, and NYGC's own GitHub repo has no fix.** Phase 3's first real benchmarking run against `COLO-829--COLO-829BL.snv.indel.final.v6.annotated.vcf` failed `tabix` indexing with parse errors. Byte-level forensic analysis (via Python's `tarfile` module, independent of the `tar` CLI, to rule out a reader bug) found the real cause: this entry's declared ~36MB of content is **not COLO-829 HiSeqX data at all** — it's the archive's other five VCF files (COLO-829 NovaSeq, both HCC-1143 platforms, both HCC1187 platforms) concatenated together, each with its own genuine leaked tar header, and the true HiSeqX content simply absent. The evidence is exact, not circumstantial: summing the other five files' real sizes (11,987,460 + 7,116,918 + 8,388,545 + 4,258,721 + 4,403,703 = 36,155,347 bytes) and adding the expected per-fragment tar-header overhead (~8,237 bytes across 4 embedded header blocks) accounts for the corrupted entry's full 36,163,584 bytes of actual content with nothing left over. A fresh re-download produced a byte-identical SHA-256 to the original, ruling out a download/transfer glitch. Independent corroboration found afterward: the GitHub repo's own README documents `Variants.HighCoverage.tar.gz` as 140MB, but the blob actually committed is only 70.2MB — almost exactly half — and a second archive in the same repo (`Variants.Downsampled.tar.gz`) shows the identical ~50% shortfall pattern. The repo has only 7 commits (all from 2026-05, evidently a from-scratch reconstruction of the original 2019 site, likely from an old web-archive snapshot), the corrupted archive has never been touched since its one commit, there are no tags/releases/alternate branches, no GitHub issue reports this (checked — 0 open, 0 closed), and the paper's own data-availability statement points back to this same GitHub-backed site with no separate mirror. **Conclusion: this specific file, as published, cannot be fixed or recovered from elsewhere — the real COLO-829 HiSeqX SNV/indel VCF is not available from NYGC's public release.** Worth reporting upstream as a GitHub issue at some point since nobody else has flagged it, but that's a courtesy to the community, not a pipeline blocker.

**Revised decision: the COLO-829 NovaSeq somatic VCF (same archive, confirmed intact) is now the primary SNV/indel truth set for Module 5 benchmarking, in place of the originally-chosen HiSeqX version.** This reopens the platform-mismatch consideration the original HiSeqX choice was specifically meant to minimize (our PRJEB27698 FASTQs are HiSeq X Ten reads, not NovaSeq — see the caveat at the end of this section) — but a real, usable truth set with one extra axis of platform difference is strictly better than zero usable data from a file that doesn't actually contain what its name claims.

The original HiSeqX-based decision record is kept below for the paper trail, followed by what actually changed for the NovaSeq switch.

**Original decision (2026-08-29, now superseded above): NYGC open somatic VCF (HiSeqX version) is the primary SNV/indel truth set for Module 5 benchmarking.**

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
  ~~Extract only the **COLO-829 (HiSeqX)** SNV/indel VCF + CNV BED, not the NovaSeq version: HiSeqX is a platform match to the PRJEB27698 FASTQs (also HiSeq X Ten), which removes one axis of difference even though sequencing centre and library prep still differ.~~ **Superseded 2026-08-30** — that HiSeqX file is corrupted at the source (see the update at the top of this section). Extract the **NovaSeq** version instead: `COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf`. Since GNU tar's own extended-header handling produced the corrupted read on the (unrelated but nearby, and now known-fine) NovaSeq entry too during earlier diagnosis, extracting via Python's `tarfile` directly is the safer, verified-working method (a plain `tar xf` may also work fine for this specific entry — it was never shown to be the problem, only the HiSeqX entry was — but `tarfile` is what was actually confirmed clean):
  ```bash
  python3 -c "
  import tarfile
  with tarfile.open('Variants.HighCoverage.tar.gz') as tf:
      member = tf.getmember('COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf')
      with tf.extractfile(member) as src, open('COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf', 'wb') as dst:
          dst.write(src.read())
  "
  ```
- **Confirmed file in hand (revised 2026-08-30):** `~/projects/somatic-variant-analysis-COLO829/truth_set/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf` — 11,987,460 bytes, uncompressed VCF (not `.vcf.gz`), verified clean (correct `##fileformat=VCFv4.2` header, single `#CHROM` line, no embedded tar artifacts — unlike the HiSeqX file above). For som.py (Module 5's benchmarking tool — not hap.py, see `docs/PHASE3_NOTES.md`) you'll want it bgzip-compressed and tabix-indexed (handled automatically by `modules/benchmarking.nf`'s `PREPARE_TRUTH_VCF`, Phase 3):
  ```bash
  bgzip -k COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf
  tabix -p vcf COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz
  ```
  (The originally-extracted, corrupted `COLO-829--COLO-829BL.snv.indel.final.v6.annotated.vcf` should be deleted or left unused — it is not valid VCF content.)
- **Companion CNV file also extracted (still HiSeqX, unaffected by this issue — confirmed clean per the full tar member listing in the Phase 3 diagnosis):** `COLO-829--COLO-829BL.cnv.annotated.v6.final.bed` (13,490 bytes) — not the truth set for Module 5, but usable later as an extra Module 7 cross-check alongside Valle-Inclan's SV truth set. A NovaSeq CNV bed (`COLO-829-NovaSeq--COLO-829BL-NovaSeq.cnv.annotated.v6.final.bed`, 16,081 bytes) also exists in the archive if platform-matching the CNV cross-check to the new SNV/indel truth-set choice is ever wanted.
- **Access:** Open, no login (GitHub-hosted, no dbGaP/EGA involved for this file)
- **Genome build:** **GRCh38**, confirmed from the paper's methods text — "Sequencing reads were aligned to the GRCh38 reference genome (1000 Genomes version)." (Exact decoy/patch level still TODO — see §2 below, must match what Mutect2/PoN/gnomAD are built against.)
- **Important clarification (resolved 2026-08-29):** the paper *separately* re-ran their pipeline on GRCh37 purely to compare against Craig et al. 2016's GRCh37-only, EGA-controlled-access truth set (EGAD00001002142, part of EGAS00001001385 — the same Craig 2016 dataset already declined below, not a new resource). That GRCh37 comparison run is not a file we use; it's a paragraph in their methods. **The GRCh38 `Variants.HighCoverage.tar.gz` file above is what we download**, and it's unaffected by the GRCh37 side-comparison. That comparison is actually a positive signal for our choice: the paper reports "we called over 98% of the Craig et al. SNVs" — i.e. NYGC's own calls were independently validated at 98% concordance against the real gold-standard Craig 2016 set, which is exactly the multi-platform consensus truth set we can't access directly (controlled-access).
- **Bundle:** `Variants.HighCoverage.tar.gz` — actual size **73,643,008 bytes (~70.2 MiB)**, confirmed via HTTP `content-length` header matching the downloaded file exactly (not a download problem in the sense of a truncated/incomplete transfer). README's stated "140MB" was originally assumed stale/wrong — **revised 2026-08-30: it's a real symptom, not stale metadata.** The repo's README describing a 140MB file when only a 70.2MB blob was ever committed is now understood as the visible fingerprint of the same ~50% data-loss event that corrupted the COLO-829 HiSeqX entry (see the update at the top of this section) — the whole archive is short by almost exactly half, not just the one entry we happened to need. Contains all three cell lines' HiSeqX + NovaSeq VCFs, CNV BEDs, and SV bedpe files — of these, only the **COLO-829 NovaSeq** SNV/indel VCF is both relevant to Module 5 and confirmed intact (the COLO-829 HiSeqX one is not).
- **High-confidence BED — resolved:** NYGC does **not** publish a callable-regions BED for SNV/indel — the archive only contains a CNV bed and SV/SV-high-confidence bedpe per cell line, no SNV/indel-specific confidence regions. Module 5's benchmarking (som.py, not hap.py — see `docs/PHASE3_NOTES.md`) will run without a confidence-region restriction unless a separate callability/mappability filter is applied — a Phase 3 decision, not a Phase 0 blocker.
- **Raw NYGC sequencing data (for reference, not used by us):** dbGaP study `phs001839.v1.p1` — https://dbgap.ncbi.nlm.nih.gov/beta/study/phs001839.v1.p1/ — controlled access, irrelevant to this pipeline since we source FASTQs from PRJEB27698 instead (§1 above)

**Caveats, carried into `docs/benchmarking_results.md`:**
1. This truth set was generated by NYGC's own pipeline from a separate sequencing run/library prep, not the exact PRJEB27698 FASTQs this pipeline aligns. It is not itself an independent multi-platform consensus (Craig et al. 2016 is that, but is EGA/dbGaP controlled-access only) — it's one lab's own somatic calls. However, per the clarification above, NYGC's own paper already cross-validated these calls at 98% concordance against Craig et al. 2016, which meaningfully strengthens confidence in using it as our benchmark. Still, frame results against it as "agreement with a GATK-family-adjacent somatic pipeline, itself independently validated against the field's gold standard" rather than as unmediated ground truth.
2. **Added 2026-08-30, platform mismatch now unavoidable:** the original choice of the HiSeqX version specifically to platform-match the PRJEB27698 FASTQs (also HiSeq X Ten) had to be abandoned because that specific file is corrupted at the source with no fix available (see above) — the truth set actually in use is the **NovaSeq** version instead, adding sequencing-platform as one more real axis of difference between our calls and the truth set, on top of the sequencing-centre/library-prep difference already noted in caveat 1. This should be stated plainly alongside any benchmarking numbers, not left implicit.

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

**`.tbi` index files:** Phase 2's `modules/reference_prep.nf` (`INDEX_VCF`) builds a tabix index for all three resource VCFs unconditionally if one doesn't already sit alongside them — see `docs/PHASE2_NOTES.md`. Whether the Broad bucket publishes `.tbi` companions alongside these files hasn't been independently confirmed from this session; it doesn't block anything either way.

**Melanoma driver-gene interval BED — RESOLVED 2026-08-30:** `data/gene_lists/dev_intervals.bed` now has real GRCh38 coordinates for all 8 genes in `melanoma_genes.tsv`, resolved via Ensembl's REST API (`https://rest.ensembl.org/lookup/symbol/homo_sapiens/<GENE>?content-type=application/json`) — each lookup's `assembly_name` confirmed `"GRCh38"` before use, not assumed. Padded ±2000bp per the file's own TODO note. Full detail and the TERT-promoter caveat: `docs/PHASE2_NOTES.md`.

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
| CNVkit | `quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0` | Plan's original choice, confirmed both by bioconda recipe research (Phase 0) and by a real `docker pull`/execution (Phase 4, 2026-08-31) — this one worked on the first try, unlike htslib and hap.py's guessed tags in Phase 3. Wired into `modules/cnvkit.nf`. |
| hap.py | `quay.io/biocontainers/hap.py:0.3.14--py27h5c5a3ab_0` (tag corrected 2026-08-30 — `0.3.15-0` doesn't exist, biocontainers tags are `<version>--<build-string>`, confirmed with a real `docker pull`) | `pkrusche/hap.py:latest` confirmed unmaintained (~9 years since last push) — switched to the actively-rebuilt biocontainers image for consistency with the rest of the pipeline's `quay.io/biocontainers/*` convention. **Note (Phase 3, 2026-08-30):** despite the image name, Module 5 actually invokes **`som.py`** from inside this same container, not `hap.py` — `som.py` ships in the identical package and is the correct tool for a GT-less somatic VCF comparison; see `docs/PHASE3_NOTES.md`. |

---

## 8. BQSR decision (carried over from the plan, unchanged)

Base Quality Score Recalibration is **deliberately not applied** before Mutect2. Per the Broad's current somatic best-practices guidance, BQSR's benefit for somatic calling is less clear than for germline, and is considered optional for WGS cell-line data without FFPE artefacts. Documented here per the plan's instruction not to leave this quietly assumed.

---

## 9. gnomAD population AF resource for VEP annotation

Same open question as item 4 (gnomAD v4.1 vs. whatever's embedded in the Broad bucket file) — Module 6's annotation step and Module 4's germline-resource filtering should cite a consistent gnomAD version; confirm and record here once both are downloaded.
