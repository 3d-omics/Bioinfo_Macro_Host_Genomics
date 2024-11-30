rule align__bcftools__call:
    input:
        crams=[MARK_DUPLICATES / f"{sample_id}.cram" for sample_id in SAMPLES],
        crais=[MARK_DUPLICATES / f"{sample_id}.cram.crai" for sample_id in SAMPLES],
        fasta=REFERENCE / f"{HOST_NAME}.fa.gz",
        fai=REFERENCE / f"{HOST_NAME}.fa.gz.fai",
    output:
        bcf=BCFTOOLS / "{region}.bcf",
    log:
        BCFTOOLS / "{region}.log",
    conda:
        "../../environments/bcftools.yml"
    resources:
        mem_mb=8 * 1024,
        runtime=24 * 60,
    shell:
        """
        ( bcftools mpileup \
            --fasta-ref {input.fasta} \
            --region {wildcards.region} \
            --output-type u \
            {input.crams} \
        | bcftools call \
            --multiallelic-caller \
            --variants-only \
            --output-type u \
        | bcftools filter \
            --include 'QUAL > 30' \
            --output-type b \
            --output {output.bcf} \
        ) 2> {log}
        """


rule align__bcftools__concat:
    input:
        [BCFTOOLS / f"{region}.bcf" for region in REGIONS],
    output:
        BCFTOOLS / f"{HOST_NAME}.vcf.gz",
    log:
        BCFTOOLS / f"{HOST_NAME}.log",
    conda:
        "../../environments/bcftools.yml"
    shell:
        """
        bcftools concat \
            --output-type z \
            --output {output} \
            {input} \
        2> {log}
        """


rule align__bcftools__all:
    input:
        BCFTOOLS / f"{HOST_NAME}.vcf.gz",
