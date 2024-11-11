include: "genotype_functions.smk"


rule variants__genotype__genotype_gvcfs:
    """Genotype a single region"""
    input:
        vcf_gz=CALL / "{region}.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
    output:
        vcf_gz=GENOTYPE / "{region}.vcf.gz",
    log:
        GENOTYPE / "{region}.log",
    conda:
        "../../environments/gatk4.yml"
    shell:
        """
        gatk GenotypeGVCFs \
            --variant {input.vcf_gz} \
            --reference {input.reference} \
            --output {output.vcf_gz} \
        2> {log} 1>&2
        """


rule variants__genotype__genotype_gvcfs__all:
    input:
        [GENOTYPE / f"{region}.vcf.gz" for region in REGIONS],


rule variants__genotype__merge_vcfs:
    """Join all the GVCFs into a single one

    Mysterioustly MergeVcfs fucks up the file
    """
    input:
        vcf_gz=[GENOTYPE / f"{region}.vcf.gz" for region in REGIONS],
    output:
        vcf_gz=GENOTYPE / "all.vcf.gz",
        tbi=GENOTYPE / "all.vcf.gz.tbi",
    log:
        GENOTYPE / "all.log",
    conda:
        "../../environments/bcftools.yml"
    params:
        input_string=compose_merge_vcfs_input_line,
    shell:
        """
        bcftools concat \
            --output {output.vcf_gz} \
            --output-type z \
            --write-index=tbi \
            {input.vcf_gz} \
        2> {log} 1>&2
        """


rule variants__genotype__merge_vcfs__all:
    input:
        vcf_gz=GENOTYPE / "all.vcf.gz",


rule variants__genotype__all:
    input:
        rules.variants__genotype__genotype_gvcfs__all.input,
        rules.variants__genotype__merge_vcfs__all.input,
