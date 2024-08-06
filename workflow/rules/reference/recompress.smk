rule reference__recompress__genome__:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=features["reference"]["dna"],
    output:
        fa_gz=REFERENCE / "genome.fa.gz",
    log:
        REFERENCE / "genome.log",
    conda:
        "__environment__.yml"
    cache: True
    shell:
        """
        ( gzip \
            --decompres \
            --stdout \
            {input.fa_gz} \
        | bgzip \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference__recompress__vcf__:
    """Extract the vcf.gz on config.yaml into known_variants.vcf.gz with bgzip"""
    input:
        vcf_gz=features["reference"]["known_vcf"],
    output:
        vcf_gz=REFERENCE / "known_variants.vcf.gz",
    log:
        REFERENCE / "known_variants.log",
    conda:
        "__environment__.yml"
    cache: True
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


rule reference__recompress__gtf__:
    """Extract the vcf.gz on config.yaml into known_variants.vcf.gz with bgzip"""
    input:
        gtf_gz=features["reference"]["gtf"],
    output:
        gtf_gz=REFERENCE / "annotation.gtf.gz",
    log:
        REFERENCE / "annotation.log",
    conda:
        "__environment__.yml"
    cache: True
    shell:
        """
        ( bedtools sort \
            -i {input.gtf_gz} \
        | bgzip \
            --threads {threads} \
        > {output.gtf_gz} \
        ) 2> {log}
        """


rule reference__recompress:
    input:
        rules.reference__recompress__genome.output,
        rules.reference__recompress__vcf.output,
        rules.reference__recompress__gtf.output,
