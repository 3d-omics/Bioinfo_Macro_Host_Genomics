rule helpers__samtools__gzfai__:
    """Index fa.gz"""
    input:
        "{prefix}.fa.gz",
    output:
        "{prefix}.fa.gz.fai",
    log:
        "{prefix}.fa.gz.fai.log",
    conda:
        "../../environments/samtools.yml"
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule helpers__samtools__crai__:
    """Generate a cram index"""
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    log:
        "{prefix}.cram.crai.log",
    conda:
        "../../environments/samtools.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule helpers__samtools__dict_fa__:
    """Generate a dictionary from a .fa"""
    input:
        "{prefix}.fa",
    output:
        "{prefix}.dict",
    log:
        "{prefix}.dict.log",
    conda:
        "../../environments/samtools.yml"
    shell:
        "samtools dict {input} --output {output} 2> {log} 1>&2"


rule helpers__samtools__dict_fagz__:
    """Generate a dictionary from a fa.gz"""
    input:
        "{prefix}.fa.gz",
    output:
        "{prefix}.dict",
    log:
        "{prefix}.dict.log",
    conda:
        "../../environments/samtools.yml"
    shell:
        "samtools dict {input} --output {output} 2> {log} 1>&2"


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


rule helpers__samtools__samtools_stats_cram__:
    """Compute stats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.stats.tsv",
    log:
        "{prefix}.stats.log",
    conda:
        "../../environments/samtools.yml"
    shell:
        "samtools stats {input.cram} > {output.tsv} 2> {log}"


rule helpers__samtools__samtools_flagstats_cram__:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "../../environments/samtools.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule helpers__samtools__samtools_idxstats_cram__:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "../../environments/samtools.yml"
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"


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
