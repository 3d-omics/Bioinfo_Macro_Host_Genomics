rule align__bcftools__call:
    input:
        crams=[MARK_DUPLICATES / f"{sample_id}.cram" for sample_id in SAMPLES],
        crais=[MARK_DUPLICATES / f"{sample_id}.cram.crai" for sample_id in SAMPLES],
        fasta=REFERENCE / "genome.fa.gz",
        fai=REFERENCE / "genome.fa.gz.fai",
    output:
        bcf=BCFTOOLS / "{region}.bcf",
    log:
        BCFTOOLS / "{region}.log",
    conda:
        "../../environments/bcftools.yml"
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
        vcfs=[BCFTOOLS / f"{region}.bcf" for region in REGIONS],
    output:
        vcf=BCFTOOLS / "known_variants.vcf.gz",
    log:
        BCFTOOLS / "known_variants.log",
    conda:
        "../../environments/bcftools.yml"
    shell:
        """
        bcftools concat \
            --output-type z \
            --output {output.vcf} \
            {input.vcfs} \
        2> {log}
        """


rule align__bcftools__all:
    input:
        vcf=BCFTOOLS / "known_variants.vcf.gz",
