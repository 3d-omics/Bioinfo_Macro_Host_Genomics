include: "mark_duplicates_functions.smk"


rule align__mark_duplicates:
    """Mark duplicates for all contigs and merging samples from different libraries"""
    input:
        cram=get_crams_for_mark_duplicates,
        reference=REFERENCE / "genome.fa.gz",
    output:
        bam=pipe(MARK_DUPLICATES / "{sample_id}.bam"),
        metrics=MARK_DUPLICATES / "{sample_id}.metrics.tsv",
    log:
        MARK_DUPLICATES / "{sample_id}.bam.log",
    conda:
        "../../environments/gatk4.yml"
    params:
        input_cram=compose_input_line_for_mark_duplicates,
    threads: 0  # Pipe! The bottleneck is in samtools compression
    group:
        "align_{sample_id}"
    shell:
        """
        mkdir --parents {output.bam}.tmp 2>> {log} 1>&2

        gatk MarkDuplicates \
            {params.input_cram} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --ASSUME_SORT_ORDER coordinate \
            --COMPRESSION_LEVEL 0 \
            --REFERENCE_SEQUENCE {input.reference} \
            --TMP_DIR {output.bam}.tmp \
        2>> {log} 1>&2

        rm -rfv {output.bam}.tmp 2>> {log} 1>&2
        """


rule align__mark_duplicates__bam_to_cram:
    """Conver MarkDuplicates from BAM to CRAM

    Note: As of 2024-11-11 it is not possible to pipe directly to samtools because of a
    rogue character sent to samtools.
    """
    input:
        bam=MARK_DUPLICATES / "{sample_id}.bam",
        reference=REFERENCE / "genome.fa.gz",
    output:
        MARK_DUPLICATES / "{sample_id}.cram",
    log:
        MARK_DUPLICATES / "{sample_id}.cram.log",
    conda:
        "../../environments/samtools.yml"
    group:
        "align_{sample_id}"
    shell:
        """
        samtools view \
            --threads {threads} \
            --output-fmt CRAM \
            --reference {input.reference} \
            --output {output} \
            {input.bam} \
        2> {log} 1>&2
        """


rule align__mark_duplicates__all:
    """Mark duplicates in all chromosomes and all libraries"""
    input:
        [MARK_DUPLICATES / f"{sample_id}.cram" for sample_id in SAMPLES],
