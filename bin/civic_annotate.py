#!/usr/bin/env python3
"""civic_annotate.py -- Module 8 (oncogenicity/actionability interpretation), added 2026-09-02.

Reads a SnpEff-annotated VCF (ANN= field per record, see modules/oncogenicity.nf) and, for each
PASS record, tries to identify a CIViC (civicdb.org) variant matching the called protein change,
then reports whatever clinical evidence CIViC has for it.

Why CIViC and not COSMIC's Cancer Gene Census: COSMIC CGC registration has been on hold since
Phase 0 (docs/PHASE0_FINDINGS.md) -- no account category cleanly fits a personal, non-commercial
portfolio project. CIViC's GraphQL API needs no registration for the query volume this module
uses (a handful of genes), confirmed live (not from docs alone) 2026-09-02:
    curl -s -X POST https://civicdb.org/api/graphql -H "Content-Type: application/json" \\
      -d '{"query":"query { gene(entrezSymbol: \\"BRAF\\") { name variants { edges { node {
      name molecularProfiles { edges { node { evidenceItems { edges { node { evidenceLevel
      evidenceType significance disease{name} therapies{name} } } } } } } } } } } }"}'
returned real data (BRAF A598V's evidence item, evidenceLevel C, disease Melanoma, therapies
listed) -- the exact field names below (`evidenceLevel`, `evidenceType`, `significance`,
`disease{name}`, `therapies{name}`) are the real, live-confirmed schema, not a guess from CIViC's
own documentation (which explicitly warns evidence items hang off molecular profiles, not
variants directly -- true, and reflected in the query shape here).

Matching strategy: CIViC names variants by protein change (e.g. "V600E", "P403fs"), not genomic
coordinate. SnpEff's ANN field gives HGVS.p in three-letter form (e.g. "p.Val600Glu"). This script
converts three-letter HGVS.p to CIViC's one-letter short form and matches by exact (case-
insensitive) name within the called gene. A record can carry multiple transcript annotations
(SnpEff reports one ANN entry per overlapping transcript); every transcript's HGVS.p is tried, in
the order SnpEff emitted them, and the first CIViC match wins -- this matters in practice for
COLO829's own BRAF V600E call, which SnpEff reports as "p.Val640Glu" on one transcript
(NM_001374258.1, an alternate numbering) and "p.Val600Glu" on another (NM_004333.6, the canonical
RefSeq transcript) -- only the second form matches CIViC's "V600E" naming.

A record with no HGVS.p at all (intronic/UTR/downstream/upstream -- MODIFIER impact, no protein
change) cannot be protein-level matched by construction, not by a matching failure -- reported
honestly as "non-coding annotation, no protein-level match applicable" rather than silently
skipped or forced into a lower-confidence gene-level guess.

Usage: civic_annotate.py <snpeff_annotated.vcf> <output.tsv>
"""
import sys
import json
import time
import re
import urllib.request
import urllib.error

CIVIC_API = "https://civicdb.org/api/graphql"

# Standard 3-letter -> 1-letter amino acid code map, plus the stop-codon/termination codes SnpEff
# and HGVS both use.
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}

# Matches HGVS.p forms this project's real Mutect2/SnpEff output actually produced:
#   p.Val600Glu   (missense: 3-letter, digits, 3-letter)
#   p.Ala68fs     (frameshift: 3-letter, digits, literal "fs")
# Anything else (in-frame del/ins, stop-gain written differently, etc.) is deliberately left
# unmatched rather than guessed at -- not encountered in this project's real calls so far, and
# CIViC's own naming for those event types varies enough that a wrong guess is worse than an
# honest "not attempted."
HGVS_P_RE = re.compile(r"^p\.([A-Za-z]{3})(\d+)(fs|[A-Za-z]{3})$")


def hgvs_p_to_civic_name(hgvs_p):
    """Convert e.g. 'p.Val600Glu' -> 'V600E', 'p.Ala68fs' -> 'A68fs'. Returns None if the form
    isn't one of the two patterns this project has actually seen (see HGVS_P_RE comment)."""
    if not hgvs_p:
        return None
    m = HGVS_P_RE.match(hgvs_p)
    if not m:
        return None
    aa1, pos, tail = m.groups()
    aa1_1letter = AA3_TO_1.get(aa1)
    if aa1_1letter is None:
        return None
    if tail == "fs":
        return f"{aa1_1letter}{pos}fs"
    aa2_1letter = AA3_TO_1.get(tail)
    if aa2_1letter is None:
        return None
    return f"{aa1_1letter}{pos}{aa2_1letter}"


def civic_query_gene(gene_symbol, max_retries=3):
    """Fetch every variant CIViC has for a gene, each with its evidence items. One query per
    gene (not per variant) -- cheap at this project's scale (a handful of genes total), and
    avoids needing to guess whether CIViC's API supports server-side filtering by variant name
    (not confirmed live, so not relied on)."""
    query = """
    query($symbol: String) {
      gene(entrezSymbol: $symbol) {
        name
        variants {
          edges {
            node {
              name
              molecularProfiles {
                edges {
                  node {
                    name
                    evidenceItems {
                      edges {
                        node {
                          evidenceLevel
                          evidenceType
                          significance
                          disease { name }
                          therapies { name }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    """
    payload = json.dumps({"query": query, "variables": {"symbol": gene_symbol}}).encode("utf-8")
    req = urllib.request.Request(
        CIVIC_API, data=payload, headers={"Content-Type": "application/json"}, method="POST"
    )
    last_err = None
    for attempt in range(1, max_retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last_err = e
            print(f"WARN: CIViC query for {gene_symbol} failed (attempt {attempt}/{max_retries}): {e}", file=sys.stderr)
            time.sleep(2 * attempt)
    print(f"FAILED: CIViC query for {gene_symbol} after {max_retries} attempts: {last_err}", file=sys.stderr)
    return None


def build_variant_index(civic_response):
    """gene response -> {variant_name.lower(): [evidence_item dict, ...]}"""
    index = {}
    if not civic_response or not civic_response.get("data") or not civic_response["data"].get("gene"):
        return index
    for v_edge in civic_response["data"]["gene"]["variants"]["edges"]:
        v = v_edge["node"]
        items = []
        for mp_edge in v["molecularProfiles"]["edges"]:
            for ei_edge in mp_edge["node"]["evidenceItems"]["edges"]:
                items.append(ei_edge["node"])
        index[v["name"].lower()] = items
    return index


def parse_ann_field(info_field):
    """Extract all ANN= subentries from a VCF INFO field as a list of 16-element lists (SnpEff's
    documented ANN field order: Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|
    Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA_pos|CDS_pos|AA_pos|Distance|Errors)."""
    for field in info_field.split(";"):
        if field.startswith("ANN="):
            raw = field[len("ANN="):]
            entries = []
            for entry in raw.split(","):
                parts = entry.split("|")
                if len(parts) >= 11:
                    entries.append(parts)
            return entries
    return []


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <snpeff_annotated.vcf> <output.tsv>", file=sys.stderr)
        sys.exit(1)

    vcf_path, out_path = sys.argv[1], sys.argv[2]

    records = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt, _qual, filt, info = fields[:8]
            if filt != "PASS":
                continue
            ann_entries = parse_ann_field(info)
            records.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "ann": ann_entries,
            })

    print(f"Loaded {len(records)} PASS record(s) from {vcf_path}", file=sys.stderr)

    # Collect the set of genes we'll actually need to query, from the highest-impact ANN entry's
    # Gene_Name per record (SnpEff already orders ANN entries by impact, most severe first, per
    # its own documented sorting -- confirmed by this project's own real output: every record's
    # first ANN entry was its most-impactful available transcript annotation).
    genes_needed = set()
    for rec in records:
        if rec["ann"]:
            genes_needed.add(rec["ann"][0][3])  # Gene_Name

    print(f"Genes to query in CIViC: {sorted(genes_needed)}", file=sys.stderr)

    civic_index_by_gene = {}
    for gene in sorted(genes_needed):
        print(f"Querying CIViC for {gene}...", file=sys.stderr)
        resp = civic_query_gene(gene)
        civic_index_by_gene[gene] = build_variant_index(resp)
        time.sleep(0.5)  # courteous pacing -- well under CIViC's anonymous rate limit at this volume

    rows = []
    header = [
        "chrom", "pos", "ref", "alt", "gene", "impact", "hgvs_c", "hgvs_p",
        "civic_variant_matched", "evidence_level", "evidence_type", "significance",
        "disease", "therapies", "note",
    ]

    for rec in records:
        if not rec["ann"]:
            rows.append([rec["chrom"], rec["pos"], rec["ref"], rec["alt"], "", "", "", "",
                         "", "", "", "", "", "", "No ANN entry from SnpEff for this record"])
            continue

        gene = rec["ann"][0][3]
        impact = rec["ann"][0][2]
        variant_index = civic_index_by_gene.get(gene, {})

        matched = False
        for ann in rec["ann"]:
            hgvs_c, hgvs_p = ann[9], ann[10]
            civic_name = hgvs_p_to_civic_name(hgvs_p)
            if not civic_name:
                continue
            evidence_items = variant_index.get(civic_name.lower())
            if evidence_items is None:
                continue
            matched = True
            if not evidence_items:
                rows.append([rec["chrom"], rec["pos"], rec["ref"], rec["alt"], gene, impact,
                             hgvs_c, hgvs_p, civic_name, "", "", "", "", "",
                             "CIViC has this variant but no evidence items curated for it yet"])
            else:
                for ei in evidence_items:
                    disease = ei["disease"]["name"] if ei.get("disease") else ""
                    therapies = "|".join(t["name"] for t in ei.get("therapies", []))
                    rows.append([rec["chrom"], rec["pos"], rec["ref"], rec["alt"], gene, impact,
                                 hgvs_c, hgvs_p, civic_name, ei.get("evidenceLevel", ""),
                                 ei.get("evidenceType", ""), ei.get("significance", ""),
                                 disease, therapies, ""])
            break  # first matching transcript's annotation wins, per module docstring

        if not matched:
            # Distinguish "genuinely non-coding, protein-level matching doesn't apply" from
            # "has a protein change but CIViC doesn't have this exact variant" -- these are very
            # different findings and collapsing them would misrepresent what was actually tried.
            any_hgvs_p = any(ann[10] for ann in rec["ann"])
            if not any_hgvs_p:
                note = f"Non-coding annotation ({impact} impact) -- no protein change to match against CIViC"
            else:
                tried = sorted(set(hgvs_p_to_civic_name(ann[10]) or ann[10] for ann in rec["ann"] if ann[10]))
                note = f"No matching CIViC variant found in {gene} for: {', '.join(tried)}"
            rows.append([rec["chrom"], rec["pos"], rec["ref"], rec["alt"], gene, impact,
                         rec["ann"][0][9], rec["ann"][0][10], "", "", "", "", "", "", note])

    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(str(x) for x in row) + "\n")

    print(f"Wrote {len(rows)} row(s) to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
