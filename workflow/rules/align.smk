include: "align/reads.smk"
include: "align/bwamem2.smk"
include: "align/mark_duplicates.smk"
include: "align/bcftools.smk"
include: "align/recalibrate.smk"
include: "align/multiqc.smk"


rule align__all:
    """Run all picard steps and get all reports"""
    input:
        rules.align__reads__all.input,
        rules.align__bwamem2__all.input,
        rules.align__mark_duplicates__all.input,
        rules.align__recalibrate__all.input,
        rules.align__multiqc__all.input,
