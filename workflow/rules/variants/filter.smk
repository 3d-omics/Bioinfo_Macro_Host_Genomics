rule variants__filter__select_variants:
    """Select only snp/indes from VCF"""
    input:
        vcf=GENOTYPE / "all.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
        tbi=GENOTYPE / "all.vcf.gz.tbi",
    output:
        vcf=FILTER / "{variant_type}.raw.vcf.gz",
    log:
        FILTER / "{variant_type}.raw.log",
    conda:
        "../../environments/gatk4.yml"
    params:
        variant_type=lambda w: w.variant_type,
    shell:
        """
        gatk SelectVariants \
            --reference {input.reference} \
            --variant {input.vcf} \
            --output {output.vcf} \
            --select-type-to-include {params.variant_type} \
        2> {log} 1>&2
        """


rule variants__filter__select_variants__all:
    input:
        [FILTER / f"{variant_type}.raw.vcf.gz" for variant_type in ["SNP", "INDEL"]],


rule variants__filter__variant_filtration:
    """Filter variants for a single chromosome"""
    input:
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        vcf=FILTER / "{variant_type}.raw.vcf.gz",
    output:
        vcf=FILTER / "{variant_type}.filtered.vcf.gz",
    log:
        FILTER / "{variant_type}.log",
    conda:
        "../../environments/gatk4.yml"
    params:
        filter_name=lambda w: w.variant_type,
        filter_expression=lambda w: params["variants"]["filter"][w.variant_type],
    shell:
        """
        gatk VariantFiltration \
            --reference {input.reference} \
            --variant {input.vcf} \
            --output {output.vcf} \
            --filter-expression '{params.filter_expression}' \
            --filter-name '{params.filter_name}' \
        2> {log} 1>&2
        """


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
