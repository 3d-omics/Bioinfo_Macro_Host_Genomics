rule reference__recompress__genome:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=features["reference"]["dna"],
    output:
        fa_gz=REFERENCE / f"{HOST_NAME}.fa.gz",
    log:
        REFERENCE / f"{HOST_NAME}.fa.log",
    conda:
        "../../environments/reference.yml"
    cache: "omit-software"
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.fa_gz} \
        | bgzip \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference__recompress__vcf:
    """Extract the vcf.gz on config.yaml into known_variants.vcf.gz with bgzip"""
    input:
        vcf_gz=features["reference"]["known_vcf"],
    output:
        vcf_gz=REFERENCE / f"{HOST_NAME}.vcf.gz",
    log:
        REFERENCE / f"{HOST_NAME}.vcf.log",
    conda:
        "../../environments/reference.yml"
    cache: "omit-software"
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.vcf_gz} \
        | bgzip \
            --threads {threads} \
        > {output.vcf_gz}) \
        2> {log}
        """


rule reference__recompress__gtf:
    """Extract the vcf.gz on config.yaml into known_variants.vcf.gz with bgzip"""
    input:
        gtf_gz=features["reference"]["gtf"],
    output:
        gtf_gz=REFERENCE / f"{HOST_NAME}.gtf.gz",
    log:
        REFERENCE / f"{HOST_NAME}.gtf.log",
    conda:
        "../../environments/reference.yml"
    cache: "omit-software"
    shell:
        """
        ( bedtools sort \
            -i {input.gtf_gz} \
        | bgzip \
            --threads {threads} \
        > {output.gtf_gz} \
        ) 2> {log}
        """


rule reference__recompress__all:
    input:
        rules.reference__recompress__genome.output,
        rules.reference__recompress__vcf.output,
        rules.reference__recompress__gtf.output,
