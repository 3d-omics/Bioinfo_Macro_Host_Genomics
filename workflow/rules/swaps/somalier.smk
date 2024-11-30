rule swaps__somalier__find_sites:
    input:
        vcf=VARIANTS / "filter" / "all.filtered.vcf.gz",
        tbi=VARIANTS / "filter" / "all.filtered.vcf.gz.tbi",
    output:
        vcf=SOMALIER / "sites.vcf.gz",
    log:
        SOMALIER / "sites.log",
    conda:
        "../../environments/somalier.yml"
    params:
        min_allele_number=5,
        min_allele_frequency=0.15,
    resources:
        mem_mb=16 * 1024,
    shell:
        """
        somalier find-sites \
            --min-AN {params.min_allele_number} \
            --min-AF {params.min_allele_frequency} \
            --output-vcf {output.vcf} \
            {input.vcf} \
        2> {log} 1>&2
        """


rule swaps__somalier__extract:
    input:
        sites=SOMALIER / "sites.vcf.gz",
        variants=VARIANTS / "filter" / "all.filtered.vcf.gz",
        reference=REFERENCE / f"{HOST_NAME}.fa.gz",
        fai=REFERENCE / f"{HOST_NAME}.fa.gz.fai",
    output:
        [SOMALIER / "extracted" / f"{sample_id}.somalier" for sample_id in SAMPLES],
    log:
        SOMALIER / "extracted.log",
    conda:
        "../../environments/somalier.yml"
    params:
        out_dir=SOMALIER / "extracted",
    shell:
        """
        somalier extract \
            --out-dir {params.out_dir} \
            --sites   {input.sites} \
            --fasta   {input.reference} \
            {input.variants} \
        2> {log} 1>&2
        """


rule swaps__somalier__relate:
    input:
        extracted=[
            SOMALIER / "extracted" / f"{sample_id}.somalier" for sample_id in SAMPLES
        ],
    output:
        html=SOMALIER / "relate.html",
        pairs=SOMALIER / "relate.pairs.tsv",
        samples=SOMALIER / "relate.samples.tsv",
    log:
        SOMALIER / "relate.log",
    conda:
        "../../environments/somalier.yml"
    params:
        output_prefix=SOMALIER / "relate",
    shell:
        """
        somalier relate \
            --output-prefix {params.output_prefix} \
            --infer \
            {input.extracted} \
        2> {log} 1>&2
        """


rule swaps__somalier__all:
    input:
        rules.swaps__somalier__find_sites.output,
        rules.swaps__somalier__extract.output,
        rules.swaps__somalier__relate.output,
