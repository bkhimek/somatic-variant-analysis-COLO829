// Reference/resource preparation for GATK steps (Module 3-4, Phase 2).
//
// Mutect2, GetPileupSummaries, and FilterMutectCalls all need a FASTA index (.fai) and a
// sequence dictionary (.dict) alongside the reference -- flagged as an open TODO since Phase 0
// (docs/PHASE0_FINDINGS.md §7 action 8) and never actually built. They also need a tabix
// index (.tbi) alongside every resource VCF (panel-of-normals, germline resource, common
// biallelic sites). None of this is remotely as expensive as the bwa-mem2 index (seconds to
// low minutes, not the ~87-110GB-RAM/22-minute AWS ordeal in docs/PHASE1_NOTES.md), but the
// same "auto-detect, else build" pattern from workflows/somatic.nf's BWA_MEM2_INDEX handling is
// used here too, in case you already have some of these files (e.g. pulled `.fai`/`.dict`
// straight from the Broad FTP per docs/data_sources.md §2, or the GATK bundle VCFs already
// ship with `.tbi` companions).

process INDEX_FASTA {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path(reference_fasta)

    output:
    path("${reference_fasta}.fai"), emit: fai

    script:
    """
    samtools faidx ${reference_fasta}
    """
}

process CREATE_SEQUENCE_DICTIONARY {
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path(reference_fasta)

    output:
    path("*.dict"), emit: dict

    script:
    // GATK's default dict-naming (input basename, extension replaced with .dict) is fragile to
    // rely on implicitly -- name it explicitly via -O instead.
    def dict_name = reference_fasta.name.replaceAll(/\.(fa|fasta|fna)$/, '') + '.dict'
    // --java-options "-Xmx6g" -- applied proactively, same reasoning as modules/mutect2.nf's
    // MUTECT2 process: gatk's default JVM heap-sizing under Docker doesn't reliably scale to the
    // container's actual memory, so every GATK invocation in this pipeline sets this explicitly.
    """
    gatk --java-options "-Xmx6g" CreateSequenceDictionary -R ${reference_fasta} -O ${dict_name}
    """
}

// Generic VCF tabix-indexer, reused for whichever of panel_of_normals / germline_resource /
// common_biallelic_sites doesn't already have a .tbi sitting next to it. tag{} uses the
// filename (not a sample_id -- these are shared reference resources, not per-sample data).
process INDEX_VCF {
    tag { vcf.name }
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path(vcf)

    output:
    tuple path(vcf), path("${vcf}.tbi"), emit: vcf_with_index

    script:
    // --java-options "-Xmx6g" -- see CREATE_SEQUENCE_DICTIONARY above / modules/mutect2.nf for why
    // this is applied to every GATK invocation in this pipeline.
    """
    gatk --java-options "-Xmx6g" IndexFeatureFile -I ${vcf} -O ${vcf}.tbi
    """
}
