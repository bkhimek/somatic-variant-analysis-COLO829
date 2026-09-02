#!/usr/bin/env bash
# extract_real_gene_panel.sh -- robust per-gene extraction of real COLO829/COLO829BL reads from
# ENA's pre-aligned GRCh37 BAMs (see docs/data_sources.md section "Real-depth data via ENA's
# pre-aligned GRCh37 BAMs"), added 2026-09-02.
#
# Why this exists, not just a single `samtools view -L <whole BED file>` call (what was tried
# first): an overnight run of that single-shot approach silently truncated partway through --
# confirmed not by the stderr error alone, but by the per-sample read counts coming back wildly
# inconsistent with expected real depth (tumour: 456,807 reads; normal: only 4,891 -- real depth
# there is ~37X vs the tumour's ~98X, so a clean run should show roughly a third as many reads,
# not 1/93rd). The likely cause is a dropped connection sometime during an unattended multi-hour
# single HTTPS stream (laptop/WSL sleep, wifi hiccup, server-side idle timeout) -- not a flaw in
# the region-extraction approach itself. Splitting into one `samtools view` call per gene, each
# independently retried, means a single dropped connection only costs that one gene's retry, and
# `samtools quickcheck` after each fetch catches a truncated/incomplete BGZF stream before it can
# silently poison the merged result the way the first attempt did.
#
# Usage:
#   bash bin/extract_real_gene_panel.sh <sample_label> <remote_bam_url> <grch37_bed> <out_dir>
#
# Example (run both from the repo root):
#   bash bin/extract_real_gene_panel.sh COLO829   https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752450/COLO829T_dedup.realigned.bam data/gene_lists/dev_intervals_grch37.bed real_data
#   bash bin/extract_real_gene_panel.sh COLO829BL https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752449/COLO829R_dedup.realigned.bam data/gene_lists/dev_intervals_grch37.bed real_data

set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <sample_label> <remote_bam_url> <grch37_bed> <out_dir>" >&2
    exit 1
fi

SAMPLE_LABEL="$1"
REMOTE_BAM="$2"
BED_FILE="$3"
OUT_DIR="$4"
MAX_RETRIES=5
RETRY_DELAY_SECS=15

mkdir -p "${OUT_DIR}/${SAMPLE_LABEL}_gene_bams"
PER_GENE_BAMS=()

while IFS=$'\t' read -r chrom start end gene; do
    # Skip the BED file's own header/comment lines (dev_intervals_grch37.bed has several,
    # explaining the coordinate provenance) -- anything starting with "#", or a blank line.
    [[ "$chrom" =~ ^# ]] && continue
    [[ -z "$chrom" ]] && continue

    # BED is 0-based, half-open; samtools' "chrom:start-end" region syntax is 1-based inclusive
    # on both ends -- so the BED start needs +1, the BED end is already numerically correct.
    region="${chrom}:$((start + 1))-${end}"
    out_bam="${OUT_DIR}/${SAMPLE_LABEL}_gene_bams/${gene}.bam"
    PER_GENE_BAMS+=("$out_bam")

    attempt=1
    until samtools view -b "$REMOTE_BAM" "$region" -o "$out_bam" && samtools quickcheck "$out_bam"; do
        if [ "$attempt" -ge "$MAX_RETRIES" ]; then
            echo "FAILED: ${gene} (${region}) after ${MAX_RETRIES} attempts -- giving up. Re-run this script; already-completed genes will just be redone (cheap -- each is a small, targeted fetch)." >&2
            exit 1
        fi
        echo "Retry ${attempt}/${MAX_RETRIES} for ${gene} (${region}) after a failed/incomplete fetch..." >&2
        attempt=$((attempt + 1))
        sleep "$RETRY_DELAY_SECS"
    done
    echo "OK: ${gene} (${region}) -> ${out_bam} ($(samtools view -c "$out_bam") reads)"
done < "$BED_FILE"

echo "All ${#PER_GENE_BAMS[@]} gene regions fetched and quickcheck-passed for ${SAMPLE_LABEL} -- merging."

MERGED_BAM="${OUT_DIR}/${SAMPLE_LABEL}_real_genepanel.bam"
samtools merge -f "$MERGED_BAM" "${PER_GENE_BAMS[@]}"
samtools sort -n -o "${OUT_DIR}/${SAMPLE_LABEL}_real_genepanel.namesorted.bam" "$MERGED_BAM"
samtools fastq \
    -1 "${OUT_DIR}/${SAMPLE_LABEL}_real_R1.fastq.gz" \
    -2 "${OUT_DIR}/${SAMPLE_LABEL}_real_R2.fastq.gz" \
    -0 /dev/null -s /dev/null -n \
    "${OUT_DIR}/${SAMPLE_LABEL}_real_genepanel.namesorted.bam"

echo "Done: ${OUT_DIR}/${SAMPLE_LABEL}_real_R1.fastq.gz / ${OUT_DIR}/${SAMPLE_LABEL}_real_R2.fastq.gz"
