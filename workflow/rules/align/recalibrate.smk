rule align__recalibrate__baserecalibrator:
    """Compute the recalibration table for a single library and chromosome"""
    input:
        cram=MARK_DUPLICATES / "{sample_id}.cram",
        crai=MARK_DUPLICATES / "{sample_id}.cram.crai",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        known_sites=REFERENCE / "known_variants.vcf.gz",
        tbi=REFERENCE / "known_variants.vcf.gz.tbi",
    output:
        table=RECALIBRATE / "{sample_id}.bsqr.txt",
    log:
        RECALIBRATE / "{sample_id}.log",
    conda:
        "../../environments/gatk4.yml"
    group:
        "align_{sample_id}"
    shell:
        """
        gatk BaseRecalibrator \
            --input {input.cram} \
            --reference {input.reference} \
            --known-sites {input.known_sites} \
            --output {output.table} \
        2> {log} 1>&2
        """


rule align__recalibrate__applybqsr:
    """Apply the recalibration table to a single library and chromosome"""
    input:
        cram=MARK_DUPLICATES / "{sample_id}.cram",
        reference=REFERENCE / "genome.fa.gz",
        table=RECALIBRATE / "{sample_id}.bsqr.txt",
        dict_=REFERENCE / "genome.dict",
    output:
        bam=pipe(RECALIBRATE / "{sample_id}.bam"),
    log:
        RECALIBRATE / "{sample_id}.log",
    conda:
        "../../environments/gatk4.yml"
    threads: 0  # pipe!
    group:
        "align_{sample_id}"
    shell:
        """
        gatk ApplyBQSR \
            --input {input.cram} \
            --reference {input.reference} \
            --bqsr-recal-file {input.table} \
            --output {output.bam} \
        2> {log} 1>&2
        """


rule align__recalibrate__bam_to_cram:
    input:
        bam=RECALIBRATE / "{sample_id}.bam",
        reference=REFERENCE / "genome.fa.gz",
    output:
        cram=RECALIBRATE / "{sample_id}.cram",
    log:
        RECALIBRATE / "{sample_id}.cram.log",
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
            --output {output.cram} \
            {input.bam} \
        2> {log} 1>&2
        """


rule align__recalibrate__all:
    """Compute recalibration for all chromosomes and libraries"""
    input:
        [RECALIBRATE / f"{sample_id}.cram" for sample_id in SAMPLES],
