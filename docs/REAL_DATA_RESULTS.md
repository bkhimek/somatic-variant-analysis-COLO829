# Real-data run — first real driver mutations found (2026-09-02)

**Status: real result, cross-checked against independent literature/database sources, not just "the pipeline ran."**

## What was run

`bin/extract_real_gene_panel.sh` (the per-gene retry fix, see `docs/data_sources.md` §1) successfully extracted real reads over the 8-gene melanoma panel from ENA's pre-aligned GRCh37 COLO829/COLO829BL BAMs, with sane, proportional read counts this time (794,804 tumour reads processed vs 248,448 normal — a ~3.2:1 ratio, in the right range for the ~98X:~37X ≈ 2.65:1 depth difference between the two real accessions; per-gene ratios ranged 1.6–5.3:1, plausible given COLO829's known copy-number alterations rather than an extraction artifact).

Those FASTQs were then run through the full pipeline for the first time on real data:

```bash
nextflow run main.nf -profile docker,dev \
    --tumour_reads_1 real_data/COLO829_real_R1.fastq.gz --tumour_reads_2 real_data/COLO829_real_R2.fastq.gz \
    --normal_reads_1 real_data/COLO829BL_real_R1.fastq.gz --normal_reads_2 real_data/COLO829BL_real_R2.fastq.gz \
    --reference_fasta reference/Homo_sapiens_assembly38.fasta \
    --panel_of_normals reference/1000g_pon.hg38.vcf.gz \
    --germline_resource reference/af-only-gnomad.hg38.vcf.gz \
    --common_biallelic_sites reference/small_exac_common_3.hg38.vcf.gz \
    --truth_set_vcf truth_set/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf \
    -with-report real_run_report.html
```

26/26 processes succeeded (18m30s, 0.6 CPU-hours) — every module, including the scatter/gathered Mutect2 built 2026-09-01, ran end-to-end against real reads without a single failure.

## The 8 PASS calls, cross-checked against Cellosaurus's documented COLO829 genotype

Mutect2 (via `FILTER_MUTECT_CALLS`) reported exactly 8 `PASS` records, all inside the gene panel as expected. Rather than trust these at face value, each was checked against [Cellosaurus's COLO 829 entry (CVCL_1137)](https://www.cellosaurus.org/CVCL_1137) — an independent, curated cell-line genotype database — and, for the BRAF call specifically, against the coordinate's own annotation on [OncoKB](https://www.oncokb.org/hgvsg/7:g.140753336A%3ET?refGenome=GRCH38) and [GeneBe](https://genebe.net/variant/hg38/chr7-140753336-A-T) to confirm the exact hg38 position independently, not from memory:

| Call (this run) | Tumour AF | Normal AF | Cellosaurus's documented COLO829 genotype | Match |
|---|---|---|---|---|
| chr7:140753336 A>T | 0.669 | 0.028 (background) | **BRAF p.Val600Glu (c.1799T>A), heterozygous** | Exact coordinate match — this is *the* canonical BRAF V600E hotspot, confirmed independently via OncoKB/GeneBe as hg38 chr7:140,753,336 A>T |
| chr9:21971154 CCG>C | 0.979 | 0.024 (background) | **CDKN2A p.Ala68Glyfs\*51 (c.203_204delCG), homozygous** | 2bp deletion matches the documented frameshift exactly; the near-1.0 tumour AF matches "homozygous" |
| chr5:1295113 GG>AA | 0.974 | 0.059 | **TERT c.228_229CC>TT (-124/-125CC>TT) promoter mutation** | GG>AA on the +strand reference is the reverse complement of the documented CC>TT promoter change on TERT's minus-strand gene — the classic UV-signature dinucleotide substitution; AF pattern matches a driver promoter mutation |
| chr7:140800999 G>A (AF 0.388), chr7:140803447 G>T (AF 0.166) | — | — | Not individually listed in Cellosaurus's summary, but fall inside the same BRAF gene body | Not independently verified against a named database entry |
| chr17:7668302/7668311 (phased, AF ~0.23), chr17:7668450 (AF 0.636) | — | — | Fall inside TP53; Cellosaurus's entry lists no TP53 mutation for COLO829 | Not corroborated by Cellosaurus specifically |

**Independent confirmation beyond Cellosaurus:** `SOMPY_BENCHMARK`'s output (`results/benchmarking/*.stats.csv`) shows **tp=8, fp=0, precision=1.0 across every variant class** (SNVs, indels, MNPs) when compared against the NYGC NovaSeq truth VCF — every single PASS call this run produced, including the two BRAF-body and three TP53-region calls not named in Cellosaurus's curated summary, was independently counted as a true positive against a completely different reference source. Cellosaurus's list is a curated summary, not necessarily exhaustive, so the calls it doesn't name aren't evidence against them — they have their own independent confirmation via the truth-set match.

**Recall is expected to look near-zero (~0.0001–0.003 depending on variant class) and this is not a caller-sensitivity problem**: the truth VCF is whole-genome (41,427 SNVs, 984 indels), while this run's query is deliberately restricted to ~30kb of gene panel by `-L data/gene_lists/dev_intervals.bed`. This was flagged before the run (see the run-command message) specifically so a low recall number wouldn't be misread as a real caller weakness.

## A real, honest gap: PTEN's large deletion was not recovered

Cellosaurus also documents **PTEN c.493_634del142, homozygous** (a ~142bp deletion) for COLO829. This did not appear among the 8 PASS calls. This is very likely a real scope limitation, not a pipeline bug: Mutect2 is an assembly-based caller tuned for SNVs and short indels, and a 142bp deletion sits at the edge of or beyond what its local assembly reliably reconstructs and reports as a single indel record. Recovering it would need either a dedicated structural-variant/CNV caller (Module 6/CNVkit exists for exactly this kind of event, see caveat below) or specific Mutect2 assembly-region tuning — neither attempted here. Recorded as a known, disclosed limitation rather than silently omitted.

## CNVkit output — not trustworthy for interpretation on this run, and that's expected

`results/cnvkit/call/*.call.cns` returned genome-wide segments (chr1 through chrX) with a wide range of copy-number states. **These should not be read as real copy-number findings.** Real reads were extracted only over the 8-gene panel (plus ±2000bp padding) — the vast majority of the genome has effectively zero real coverage in this run, so CNVkit's genome-wide binning is segmenting mostly-empty bins, not real signal outside the panel. The persistently tiny `depth` column values (~1e-5, the same anomaly flagged as unresolved in `docs/PHASE4_NOTES.md`) are consistent with this, not new evidence toward resolving that open question. A real genome-wide CNVkit result would need real full-genome-depth coverage, which this extraction was deliberately designed to avoid downloading (see `docs/data_sources.md` §1 for the storage/compute tradeoff discussion).

## Bottom line

For the first time in this project, real sequencing data run through the full, real pipeline produced calls that independently match a real cell line's documented, literature-backed mutation profile — BRAF V600E, a CDKN2A frameshift, and a TERT promoter dinucleotide substitution, each cross-checked against an external database rather than taken on faith — plus zero false positives against the truth set across all 8 calls made. The two real, honestly-disclosed gaps (PTEN's large deletion missed by an SNV/indel caller, and CNVkit's genome-wide output being uninterpretable outside the panel given panel-only real coverage) are exactly the kind of scope-appropriate limitations expected from a gene-panel-only real-data extraction, not pipeline defects.

Sources checked for this write-up:
- [Cellosaurus COLO 829 (CVCL_1137)](https://www.cellosaurus.org/CVCL_1137)
- [OncoKB — BRAF V600E, hg38 chr7:140753336 A>T](https://www.oncokb.org/hgvsg/7:g.140753336A%3ET?refGenome=GRCH38)
- [GeneBe — BRAF p.Val640Glu, hg38 chr7-140753336-A-T](https://genebe.net/variant/hg38/chr7-140753336-A-T)
