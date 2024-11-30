rule annotate__snpeff__download:
    """Download a SNPEff database"""
    output:
        directory(SNPEFF_DB / "{snpeff_db}"),
    log:
        SNPEFF_DB / "{snpeff_db}.log",
    params:
        reference=lambda w: w.snpeff_db,
    wrapper:
        "v5.2.1/bio/snpeff/download"


rule annotate__snpeff__annotate:
    """Annotate variants with a SNPEff database"""
    input:
        calls=FILTER / "all.filtered.vcf.gz",
        db=SNPEFF_DB / "{snpeff_db}",
    output:
        calls=SNPEFF / "variants_{snpeff_db}.vcf.gz",
        # genes=SNPEFF / "snpEff_stats_{snpeff_db}.genes.txt",
        stats=SNPEFF / "snpEff_summary_{snpeff_db}.html",
        csvstats=SNPEFF / "snpEff_stats_{snpeff_db}.csv",
    log:
        SNPEFF / "variants_{snpeff_db}.log",
    wrapper:
        "v5.2.1/bio/snpeff/annotate"


rule annotate__snpeff__all:
    """Run SNPeff for all databases"""
    input:
        [SNPEFF / f"snpEff_stats_{snpeff_db}.csv" for snpeff_db in SNPEFF_DBS],
