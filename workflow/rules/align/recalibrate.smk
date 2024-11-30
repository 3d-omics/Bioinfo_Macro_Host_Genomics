rule align__recalibrate__baserecalibrator:
    """Compute the recalibration table for a single library and chromosome"""
    input:
        bam=MARK_DUPLICATES / "{sample_id}.cram",
        crai=MARK_DUPLICATES / "{sample_id}.cram.crai",
        ref=REFERENCE / f"{HOST_NAME}.fa.gz",
        dict=REFERENCE / f"{HOST_NAME}.dict",
        known=REFERENCE / f"{HOST_NAME}.vcf.gz",
        tbi=REFERENCE / f"{HOST_NAME}.vcf.gz.tbi",
    output:
        recal_table=RECALIBRATE / "{sample_id}.bsqr.txt",
    log:
        RECALIBRATE / "{sample_id}.log",
    group:
        "align_{sample_id}"
    wrapper:
        "v5.2.1/bio/gatk/baserecalibrator"


rule align__recalibrate__applybqsr:
    """Apply the recalibration table to a single library and chromosome"""
    input:
        cram=MARK_DUPLICATES / "{sample_id}.cram",
        reference=REFERENCE / f"{HOST_NAME}.fa.gz",
        table=RECALIBRATE / "{sample_id}.bsqr.txt",
        dict_=REFERENCE / f"{HOST_NAME}.dict",
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
        reference=REFERENCE / f"{HOST_NAME}.fa.gz",
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
