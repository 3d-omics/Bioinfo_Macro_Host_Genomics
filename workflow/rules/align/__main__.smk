include: "reads.smk"
include: "bwamem2.smk"
include: "mark_duplicates.smk"
include: "bcftools.smk"
include: "recalibrate.smk"
include: "multiqc.smk"


rule align:
    """Run all picard steps and get all reports"""
    input:
        rules.align__recalibrate.input,
        rules.align__multiqc__all.input,
