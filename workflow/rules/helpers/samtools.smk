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


rule helpers__samtools__tabix_gtf_gz__:
    input:
        gtf_gz="{prefix}.gtf.gz",
    output:
        tbi="{prefix}.gtf.gz.tbi",
    log:
        "{prefix}.gtf.gz.tbi.log",
    conda:
        "../../environments/samtools.yml"
    shell:
        "tabix -p gff {input.gtf_gz} 2> {log}"
