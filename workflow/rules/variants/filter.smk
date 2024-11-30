rule variants__filter__select_variants:
    """Select only snp/indes from VCF"""
    input:
        vcf=GENOTYPE / "all.vcf.gz",
        tbi=GENOTYPE / "all.vcf.gz.tbi",
        ref=REFERENCE / f"{HOST_NAME}.fa.gz",
    output:
        vcf=FILTER / "{variant_type}.raw.vcf.gz",
    log:
        FILTER / "{variant_type}.raw.log",
    params:
        extra=lambda w: f"--select-type-to-include {w.variant_type}",
    wrapper:
        "v5.2.1/bio/gatk/selectvariants"


rule variants__filter__select_variants__all:
    input:
        [FILTER / f"{variant_type}.raw.vcf.gz" for variant_type in ["SNP", "INDEL"]],


rule variants__filter__variant_filtration:
    """Filter variants for a single chromosome"""
    input:
        vcf=FILTER / "{variant_type}.raw.vcf.gz",
        ref=REFERENCE / f"{HOST_NAME}.fa.gz",
        dict_=REFERENCE / f"{HOST_NAME}.dict",
        # fai=REFERENCE / f"{HOST_NAME}.fa.gz.fai",
        gzi=REFERENCE / f"{HOST_NAME}.fa.gz.gzi",
    output:
        vcf=FILTER / "{variant_type}.filtered.vcf.gz",
    log:
        FILTER / "{variant_type}.log",
    params:
        filters=lambda w: {w.variant_type: params["variants"]["filter"][w.variant_type]},
    wrapper:
        "v5.2.1/bio/gatk/variantfiltration"


rule variants__filter__variant_filtration__all:
    input:
        [
            FILTER / f"{variant_type}.filtered.vcf.gz"
            for variant_type in ["SNP", "INDEL"]
        ],


rule variants__filter__merge_vcfs:
    """Merge all VCF chromosomes"""
    input:
        snps=FILTER / "SNP.filtered.vcf.gz",
        indels=FILTER / "INDEL.filtered.vcf.gz",
    output:
        FILTER / "all.filtered.vcf.gz",
    log:
        FILTER / "all.filtered.log",
    conda:
        "../../environments/gatk4.yml"
    shell:
        """
        gatk MergeVcfs \
            --INPUT {input.snps} \
            --INPUT {input.indels} \
            --OUTPUT {output} \
        2> {log} 1>&2
        """


rule variants__filter__merge_vcfs__all:
    input:
        FILTER / "all.filtered.vcf.gz",


rule variants__filter__all:
    input:
        rules.variants__filter__select_variants__all.input,
        rules.variants__filter__variant_filtration__all.input,
        rules.variants__filter__merge_vcfs__all.input,
