rule align__recalibrate__baserecalibrator__:
    """Compute the recalibration table for a single library and chromosome"""
    input:
        cram=MARK_DUPLICATES / "{sample}.cram",
        crai=MARK_DUPLICATES / "{sample}.cram.crai",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        known_sites=REFERENCE / "known_variants.vcf.gz",
        tbi=REFERENCE / "known_variants.vcf.gz.tbi",
    output:
        table=RECALIBRATE / "{sample}.bsqr.txt",
    log:
        RECALIBRATE / "{sample}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        gatk BaseRecalibrator \
            --input {input.cram} \
            --reference {input.reference} \
            --known-sites {input.known_sites} \
            --output {output.table} \
        2> {log} 1>&2
        """


rule align__recalibrate__applybqsr__:
    """Apply the recalibration table to a single library and chromosome"""
    input:
        cram=MARK_DUPLICATES / "{sample}.cram",
        reference=REFERENCE / "genome.fa.gz",
        table=RECALIBRATE / "{sample}.bsqr.txt",
        dict_=REFERENCE / "genome.dict",
    output:
        bam=pipe(RECALIBRATE / "{sample}.bam"),
    log:
        RECALIBRATE / "{sample}.log",
    conda:
        "__environment__.yml"
    threads: 0  # pipe!
    shell:
        """
        gatk ApplyBQSR \
            --input {input.cram} \
            --reference {input.reference} \
            --bqsr-recal-file {input.table} \
            --output {output.bam} \
        2> {log} 1>&2
        """


rule align__recalibrate__bam_to_cram__:
    input:
        bam=RECALIBRATE / "{sample}.bam",
        reference=REFERENCE / "genome.fa.gz",
    output:
        cram=RECALIBRATE / "{sample}.cram",
    log:
        RECALIBRATE / "{sample}.cram.log",
    conda:
        "__environment__.yml"
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


rule align__recalibrate:
    """Compute recalibration for all chromosomes and libraries"""
    input:
        [RECALIBRATE / f"{sample}.cram" for sample in SAMPLES],
