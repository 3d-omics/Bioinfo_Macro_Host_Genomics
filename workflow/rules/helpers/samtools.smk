rule helpers__samtools__vcf_gz__:
    """bgzip a vcf file"""
    input:
        "{prefix}.vcf",
    output:
        "{prefix}.vcf.gz",
    log:
        "{prefix}.vcf.gz.log",
    conda:
        "../../environments/samtools.yml"
    shell:
        "bgzip {input} 2> {log} 1>&2"
