# Phase 1 — QC and Alignment (Modules 1–2): implementation notes

**Date:** 2026-08-29
**Status:** Code written, not yet run. This was built in a Cowork sandbox with no Nextflow, no Docker, and no access to the real FASTQs or reference. I tried to install Nextflow here to actually lint/run it (`get.nextflow.io` and GitHub releases are both blocked by this sandbox's network policy) and couldn't, so this has only been checked by careful manual read-through, not executed. Manual review did catch two real DSL2 channel-cardinality bugs before delivery (both classic Nextflow gotchas, both fixed — see `workflows/somatic.nf` comments at the `.flatMap` and `BWA_MEM2_INDEX.out.index.first()` lines for what they were), but manual review is not the same guarantee as an actual run. Treat the first real `nextflow run` on your machine as the actual test, not this write-up — there could still be issues review alone doesn't catch (a typo in an unexercised code path, a container that doesn't actually have a tool where I assumed it would, etc).

---

## What's implemented

- `modules/fastqc.nf` — `FASTQC` (per-sample) + `MULTIQC` (aggregates tumour + normal into one report — this cross-sample aggregation is Phase 1's one explicit addition beyond reusing Project 4's QC step, per plan §11)
- `modules/alignment.nf` — `BWA_MEM2_INDEX`, `BWA_MEM2_ALIGN`, `SAMTOOLS_SORT`, `MARK_DUPLICATES`, chained per plan §6 Module 2 ("BWA-MEM2 → samtools sort → GATK MarkDuplicates")
- `workflows/somatic.nf` — wires both modules together; runs the whole chain once per sample (tumour, normal), matching the plan's explicit note that Module 2 "runs twice"
- `main.nf` — entry point; builds the two-sample channel from CLI params, fails fast with a clear message if required params are missing

Not implemented here, deliberately deferred to their own phases per the plan's build order: contamination estimation (Module 3, Phase 2), Mutect2 (Module 4, Phase 2), everything after.

---

## Things you need to check/fix before this actually runs

**1. `samtools` container tag is unconfirmed.** `modules/alignment.nf`'s `SAMTOOLS_SORT` process uses `quay.io/biocontainers/samtools:1.21--h50ea8bc_0`. I could not browse quay.io's tag list from this sandbox (blocked by robots.txt, same limitation hit for CNVkit/hap.py in Phase 0) to confirm this exact build-hash suffix is real. Before running:
```bash
docker pull quay.io/biocontainers/samtools:1.21--h50ea8bc_0
```
If that fails, check https://quay.io/repository/samtools/samtools?tab=tags in your browser for the real current tag and update the `container` line in `modules/alignment.nf`.

**2. No combined bwa-mem2+samtools container is used, by design.** I looked at how nf-core/modules does this (they use a Seqera Wave-hosted mulled container with a hash-only name) and deliberately didn't copy that pattern — it doesn't match this pipeline's existing convention of one plain `quay.io/biocontainers/<tool>:<version>` per tool everywhere else. Instead, `BWA_MEM2_ALIGN` writes a `.sam` file to disk and a separate `SAMTOOLS_SORT` process (different container) sorts it. Costs one extra intermediate file per sample; keeps every container single-purpose. Reasonable to revisit later if disk I/O for the intermediate SAM becomes a real bottleneck at full-WGS scale.

**3. The `dev` profile does not make alignment itself fast.** This is worth understanding before you run `-profile dev` expecting a quick smoke test. `params.interval_list` only restricts steps that take a GATK `-L` argument (Mutect2 and friends, arriving in Phase 2) — it has no effect on `BWA_MEM2_ALIGN`, which has no genomic coordinates to restrict on before alignment exists. If you want a fast Phase 1 smoke test, point `tumour_reads_1/2` and `normal_reads_1/2` at a small subsampled FASTQ pair yourself (e.g. `seqtk sample` a few thousand reads, or extract reads from a chr21 region of an existing BAM if you have one) — the pipeline won't do this for you. Full explanation is in the comment block above the `dev` profile in `nextflow.config`.

**4. BWA-MEM2 indexing runs inside the pipeline, not as a separate manual step.** `BWA_MEM2_INDEX` builds the index from whatever `params.reference_fasta` points at, every run (Nextflow's `-resume` cache will skip re-running it once it's succeeded once for that exact reference file). This means the first run against the full GRCh38 reference will spend real time and disk on indexing before alignment even starts — budget for that in addition to the alignment/sort/dedup time itself.

**5. Reference `.fai`/`.dict` are not generated yet.** Only the BWA-MEM2 index is built in Phase 1, since that's all Module 2 needs. GATK steps in Phase 2 (Mutect2, GetPileupSummaries) will need a `.fai` and a sequence dictionary (`.dict`) — add a `samtools faidx` + `gatk CreateSequenceDictionary` step when Phase 2 starts, don't assume it's already there.

**6. Sample naming is hardcoded in `main.nf`**: `COLO829` (tumour) / `COLO829BL` (normal). This matches the naming already used throughout `docs/` and `docs/run_manifest.example.json` — keep it consistent going into Phase 2, since Mutect2's tumour/normal pairing and the run manifest's `samples` block both key off these exact strings.

---

## How to run it (once you've checked the items above)

Minimal smoke test — point every required param at real files on your machine:

```bash
conda activate nextflow
cd ~/projects/somatic-variant-analysis-COLO829

nextflow run main.nf -profile docker,dev \
    --tumour_reads_1  /path/to/COLO829_R1.fastq.gz \
    --tumour_reads_2  /path/to/COLO829_R2.fastq.gz \
    --normal_reads_1  /path/to/COLO829BL_R1.fastq.gz \
    --normal_reads_2  /path/to/COLO829BL_R2.fastq.gz \
    --reference_fasta /path/to/GRCh38.fasta \
    -resume
```

If it fails on the `SAMTOOLS_SORT` container pull, that's item 1 above — fix the tag and rerun with `-resume` (BWA-MEM2 indexing and alignment results already computed will be reused, not redone).

Expected outputs once it completes:
- `results/qc/multiqc/multiqc_report.html` — combined tumour+normal QC report
- `results/alignment/COLO829/COLO829.dedup.bam` (+ `.bai`, + `.dedup.metrics.txt`)
- `results/alignment/COLO829BL/COLO829BL.dedup.bam` (+ `.bai`, + `.dedup.metrics.txt`)

Per the plan's Phase 1 verification step ("verify duplicate rates and coverage depth are appropriate for WGS somatic calling"): check the `.dedup.metrics.txt` `PERCENT_DUPLICATION` field and eyeball coverage via `samtools depth` or the MultiQC report before moving to Phase 2 — a pipeline that runs to completion isn't the same as one whose output looks right.
