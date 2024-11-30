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
        # tbi=GENOTYPE / "{region}.vcf.gz.tbi",
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
        calls=[GENOTYPE / f"{region}.vcf.gz" for region in REGIONS],
    output:
        vcf_gz=GENOTYPE / "all.vcf.gz",
        # tbi=GENOTYPE / "all.vcf.gz.tbi",
    log:
        GENOTYPE / "all.log",
    wrapper:
        "v5.2.1/bio/bcftools/concat"


rule variants__genotype__merge_vcfs__all:
    input:
        vcf_gz=GENOTYPE / "all.vcf.gz",


rule variants__genotype__all:
    input:
        rules.variants__genotype__genotype_gvcfs__all.input,
        rules.variants__genotype__merge_vcfs__all.input,
