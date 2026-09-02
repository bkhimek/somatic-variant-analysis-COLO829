# Module 8 — Oncogenicity/actionability interpretation (2026-09-02)

**Status: every component, including `bin/civic_annotate.py` itself, manually validated against real data; not yet run as a wired Nextflow DAG.** Every real risk this design raised was checked against a live system before being coded around, not guessed — see below. What's still open is running `nextflow run main.nf ...` itself with this module wired in; that's a genuinely new step, distinct from the manual `docker run`/standalone-script testing that validated each piece individually.

## Why CIViC, not COSMIC's Cancer Gene Census

The project plan's original Module 8 source is COSMIC's CGC. Phase 0 (`docs/PHASE0_FINDINGS.md`) found registration doesn't cleanly fit a personal, non-commercial portfolio project (no account category fits a private individual) and left it on hold, documenting CIViC (civicdb.org) as the open, no-registration fallback. That fallback is what's actually used here.

## Design fork: protein-level matching vs. gene-level context only

CIViC names variants by protein change ("V600E"), not genomic coordinate. Matching a Mutect2 call to CIViC evidence needs an intermediate step translating "chr7:140753336 A>T" into "p.Val600Glu" first. Two options were on the table, decided with the user 2026-09-02: gene-level-only CIViC lookup (simpler, one dependency, but only tells you a gene is relevant — not that this exact mutation is actionable), or full protein-level matching via SnpEff annotation (more moving parts, but the only version that's genuinely reusable on a future variant this project hasn't already manually cross-referenced). **Chose full SnpEff matching.**

## Every real risk checked before being coded around

This sandbox cannot reach `civicdb.org` or `quay.io` directly (both return 403 at its egress proxy, the same restriction pattern as the Broad FTP/`storage.googleapis.com` blocks from Phase 1/2) — every check below ran on the user's WSL machine, live, not guessed from documentation:

1. **Container tag.** `quay.io/biocontainers/snpeff` has many tags; picking the wrong one has already cost a round-trip once in this project (Phase 3's hap.py tag). Queried quay.io's tag API and sorted by actual `last_modified` date (not tag-string order, which can mislead — a numerically-higher-looking tag isn't necessarily newer) — `5.4.0c--hdfd78af_0` (04 Mar 2026) came back as the real newest build. Confirmed pullable.

2. **Database naming convention.** `Homo_sapiens_assembly38.fasta` uses `chr`-prefixed contigs (`docs/data_sources.md` §2). SnpEff ships both Ensembl-convention databases (`GRCh38.99`, `GRCh38.86`, etc. — bare contig names, "1" not "chr1") and a UCSC/RefSeq-convention `hg38` database. Using an Ensembl-named database against a `chr`-prefixed reference would hit the exact bare-vs-`chr` mismatch this project already found once with the ENA GRCh37 BAMs (`docs/data_sources.md` §1). Rather than guess, ran the real annotation against the real chr7:140753336 BRAF call and confirmed `hg38` produces a correct, non-error result: `NM_004333.6:c.1799T>A:p.Val600Glu` — an exact match to both Cellosaurus's documented COLO829 genotype and to OncoKB/ClinVar's independent annotation of that same hg38 coordinate (checked earlier, `docs/REAL_DATA_RESULTS.md`).

3. **SnpEff's own database-hosting has a real failure history.** A GitHub issue (pcingola/SnpEff#596) documents downloads failing from `snpeff.blob.core.windows.net`. The URLs `snpEff databases` actually returned in this session point at `snpeff-public.s3.amazonaws.com` instead — suggesting a migration since that issue — and the real download in this session succeeded cleanly. Not proof it'll always succeed; if `SNPEFF_DOWNLOAD` fails at runtime with a download error, that's this known upstream risk resurfacing, not a pipeline bug, and a retry (or checking SnpEff's GitHub issues for current status) is the right response, not assuming the module is broken.

4. **A real `OutOfMemoryError`, found and fixed.** The first real annotation attempt ran SnpEff directly against the full `filtered.vcf.gz` (276 candidate records — every filter-failed multiallelic/germline/contamination call Mutect2 considered, not just the 8 `PASS` ones) and crashed partway through with `java.lang.OutOfMemoryError`, after two other real parse errors on complex multiallelic sites. Since this module only ever needs the `PASS` calls, the fix — filtering to `PASS` first (`FILTER_PASS_VARIANTS`, plain `awk` on the FILTER column, since `bcftools` turned out not to be installed on the user's machine either) — is the correct scope, not a workaround for a limitation being lived with. Re-running against the 8-record PASS-only file completed cleanly, zero errors.

5. **CIViC's live GraphQL schema**, confirmed by a real query (not from CIViC's own docs, which explicitly warn evidence items hang off molecular profiles rather than variants directly — true, and reflected below):
   ```graphql
   query { gene(entrezSymbol: "BRAF") { name variants { edges { node {
     name molecularProfiles { edges { node {
       evidenceItems { edges { node {
         evidenceLevel evidenceType significance disease{name} therapies{name}
       } } }
     } } }
   } } } } }
   ```
   Returned real data — BRAF A598V's evidence item (level C, disease Melanoma, therapies BRAF Inhibitor / MEK Inhibitor) — confirming every field name used in `bin/civic_annotate.py` is real, not guessed. No registration or API key needed at this project's query volume (a handful of genes).

## `bin/civic_annotate.py` itself, tested standalone before any Nextflow run (2026-09-02)

Ran directly (`python3 bin/civic_annotate.py results/mutect2/annotated_pass.vcf report.tsv`) against the real 8-record annotated PASS VCF, before wiring into a Nextflow DAG. All 8 behaved exactly as designed:

- **BRAF V600E** matched CIViC's `V600E` exactly — 132 real evidence rows returned (evidence levels B/C/D, spanning melanoma sensitivity/resistance, thyroid cancer, colorectal cancer, and several other indications). The large row count is real breadth of curation for one of oncology's most-studied variants, not a bug.
- **TERT promoter, both BRAF-intronic calls, all three TP53 downstream/UTR calls** (5 records) — correctly reported `civic_variant_matched` empty with the honest note "Non-coding annotation (MODIFIER impact) — no protein change to match against CIViC." Exactly the intended behavior for calls that categorically have no protein-level HGVS.p.
- **CDKN2A frameshift** — tried all three real candidate names generated from its three annotated transcripts (`A17fs`, `A68fs`, `G83fs`) against CIViC's real CDKN2A variant list (confirmed live: `Q70fs*78`, `S56fs*51`, plus broad category variants like `Loss`/`Deletion`/`Mutation` this module deliberately does not fall back to, since matching to a category isn't evidence for the specific called variant) — no exact match, correctly reported as such rather than guessed. Worth noting, not chasing further right now: CIViC's `S56fs*51` shares the exact `*51` stop-position suffix with Cellosaurus's documented COLO829 genotype (`p.Ala68Glyfs*51`) despite the different amino-acid position — plausibly the same real `c.203_204delCG` deletion named via CDKN2A's alternate p14ARF reading frame (CDKN2A is a known special case: two overlapping isoforms translate the same DNA change into different protein sequences), but that's a hypothesis, not confirmed, and the module correctly declines to guess at cross-isoform name equivalence rather than risk a false match here or on some other gene later.

## What the module actually does

`modules/oncogenicity.nf`, wired into `workflows/somatic.nf` as Module 8, right after `FILTER_MUTECT_CALLS`:

1. `FILTER_PASS_VARIANTS` — `awk`+`bgzip` down to PASS-only records (samtools biocontainer, already trusted elsewhere in this pipeline for bgzip/tabix).
2. `SNPEFF_DOWNLOAD` (auto-skipped if `params.snpeff_data_dir/hg38` already exists, same pattern as the `.fai`/`.dict` auto-detect) — pulls SnpEff's `hg38` database once, reused across runs.
3. `SNPEFF_ANNOTATE` — annotates the PASS-only VCF, adding the `ANN=` field (gene, impact, HGVS.c, HGVS.p per overlapping transcript).
4. `CIVIC_ANNOTATE` (`bin/civic_annotate.py`, `python:3.11-slim` — an official Python image rather than a `quay.io/biocontainers/*` one, deliberately: this step is plain HTTP+JSON, no bioinformatics tool involved) — for each PASS record, converts each transcript's HGVS.p from three-letter (`p.Val600Glu`) to CIViC's one-letter short form (`V600E`), queries CIViC once per gene needed, and reports whatever evidence items match. Output: `results/oncogenicity/civic_report.tsv`.

**Multi-transcript handling, real and not hypothetical:** COLO829's own BRAF V600E call is reported by SnpEff as `p.Val640Glu` on one transcript (`NM_001374258.1`, alternate numbering) and `p.Val600Glu` on another (`NM_004333.6`, the canonical RefSeq transcript) — only the second form matches CIViC's `V600E` naming. `bin/civic_annotate.py` tries every transcript's HGVS.p in the order SnpEff emits them (already impact-sorted, most severe first, confirmed from this project's own real output) and takes the first CIViC match, rather than trusting only the first-listed transcript.

**Disclosed scope limit:** the HGVS.p parser only recognizes two shapes — missense (`p.Val600Glu`) and frameshift (`p.Ala68fs`) — because those are the only two this project's real Mutect2 output has actually produced. A future in-frame indel or stop-gain call falls through to "no protein-level match attempted" rather than a silently wrong guess. Confirmed against this project's own real 8 PASS calls (manually, via `docker run` — see below) that the CDKN2A frameshift call (`p.Ala68fs`, matching Cellosaurus's documented `c.203_204delCG` exactly) and the TERT promoter call correctly annotate with real HGVS.c/p, and that the 5 non-coding calls (2 BRAF intronic, 3 TP53 downstream/UTR) correctly report "non-coding annotation, no protein-level match applicable" rather than a forced guess.

**Correction to `docs/REAL_DATA_RESULTS.md`'s earlier characterization:** that document described the two extra BRAF-region and three TP53-region PASS calls only as "falling inside" those genes, without qualifying impact. SnpEff's real annotation shows all five are `MODIFIER`-impact (intron_variant / downstream_gene_variant / 3_prime_UTR_variant) — non-coding, not additional missense/coding mutations. This doesn't change their status as real, truth-set-confirmed true positives (som.py's tp=8/fp=0 stands), just their functional characterization — worth knowing when interpreting them, since a non-coding passenger call and a coding driver call carry very different weight.

## Real components validated manually, before wiring into Nextflow

Every step above was run directly via `docker run` against this project's own real data on the user's WSL machine, in this session, before being wired into `modules/oncogenicity.nf` — this is what "manually validated component-by-component" means concretely:
- `snpEff databases` (real container, real tag) — confirmed `hg38` exists and confirmed the historical database-hosting risk's current URLs.
- `snpEff download -dataDir ... hg38` — real download, succeeded.
- `snpEff -dataDir ... hg38 pass_only.vcf.gz` — real annotation of the real 8 PASS calls, zero errors, correct HGVS.c/p output confirmed against Cellosaurus/OncoKB.
- The CIViC GraphQL query itself — real HTTP request, real response, real field names confirmed.

**Not yet done:** running `bin/civic_annotate.py` itself (the actual matching/reporting logic) against real data, and running the whole thing as an actual Nextflow DAG (`FILTER_PASS_VARIANTS` → `SNPEFF_DOWNLOAD`/`SNPEFF_ANNOTATE` → `CIVIC_ANNOTATE`, with the auto-detect branching in `workflows/somatic.nf`). Recommended next step: run `bin/civic_annotate.py` standalone first (cheap, no Docker/Nextflow needed — plain `python3`, against the `results/mutect2/annotated_pass.vcf` this session already produced by hand) to catch any real bug in the matching logic before spending a full pipeline run on it.
