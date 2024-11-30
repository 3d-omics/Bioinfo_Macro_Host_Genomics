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
        bam=MARK_DUPLICATES / "{sample_id}.cram",
        ref=REFERENCE / f"{HOST_NAME}.fa.gz",
        dict=REFERENCE / f"{HOST_NAME}.dict",
        recal_table=RECALIBRATE / "{sample_id}.bsqr.txt",
    output:
        bam=RECALIBRATE / "{sample_id}.cram",
    log:
        RECALIBRATE / "{sample_id}.log",
    group:
        "align_{sample_id}"
    wrapper:
        "v5.2.1/bio/gatk/applybqsr"


rule align__recalibrate__all:
    """Compute recalibration for all chromosomes and libraries"""
    input:
        [RECALIBRATE / f"{sample_id}.cram" for sample_id in SAMPLES],
