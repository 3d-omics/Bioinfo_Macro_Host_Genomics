rule annotate__snpeff__download:
    """Download a SNPEff database"""
    output:
        ddirectory(SNPEFF_DB / "{snpeff_db}"),
    log:
        SNPEFF_DB / "{snpeff_db}.log",
    params:
        reference=lambda w: w.snpeff_db,
    wrapper:
        "v5.2.1/bio/snpeff/download"


rule annotate__snpeff__annotate:
    """Annotate variants with a SNPEff database"""
    input:
        vcf=FILTER / "all.filtered.vcf.gz",
        db=SNPEFF_DB / "{snpeff_db}",
    output:
        vcf=SNPEFF / "variants_{snpeff_db}.vcf.gz",
        genes=SNPEFF / "snpEff_stats_{snpeff_db}.genes.txt",
        html=SNPEFF / "snpEff_summary_{snpeff_db}.html",
        csv=SNPEFF / "snpEff_stats_{snpeff_db}.csv",
    log:
        SNPEFF / "variants_{snpeff_db}.log",
    conda:
        "__environment__.yml"
    params:
        snpeff_db="{snpeff_db}",
        datadir=SNPEFF_DB,
        html="snpEff_summary.html",
    shell:
        """
        ( snpEff ann \
            {params.snpeff_db} \
            -dataDir {params.datadir} \
            -csvStats {output.csv} \
            -verbose \
            -i vcf \
            -o gatk \
            {input.vcf} \
        | bgzip \
            --compress-level 9 \
            --stdout \
        > {output.vcf} \
        ) 2> {log} 1>&2

        mv {params.html} {output.html} 2>> {log} 1>&2
        """


rule annotate__snpeff__all:
    """Run SNPeff for all databases"""
    input:
        [SNPEFF / f"snpEff_stats_{snpeff_db}.csv" for snpeff_db in SNPEFF_DBS],
