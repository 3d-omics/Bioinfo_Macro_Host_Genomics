include: "genotype_functions.smk"


rule variants__genotype__genotype_gvcfs:
    """Genotype a single region"""
    input:
        gvcf=CALL / "{region}.vcf.gz",
        ref=REFERENCE / f"{HOST_NAME}.fa.gz",
        # dict_=REFERENCE / f"{HOST_NAME}.dict",
        # fai=REFERENCE / f"{HOST_NAME}.fa.gz.fai",
        # gzi=REFERENCE / f"{HOST_NAME}.fa.gz.gzi",
    output:
        vcf=GENOTYPE / "{region}.vcf.gz",
        tbi=GENOTYPE / "{region}.vcf.gz.tbi",
    log:
        GENOTYPE / "{region}.log",
    wrapper:
        "v5.2.1/bio/gatk/genotypegvcfs"


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
