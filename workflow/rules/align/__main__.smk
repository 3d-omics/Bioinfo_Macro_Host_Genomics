include: "__functions__.smk"
include: "reads.smk"
include: "index.smk"
include: "map.smk"
include: "mark_duplicates.smk"
include: "bcftools.smk"
include: "recalibrate.smk"


rule align:
    """Run all picard steps and get all reports"""
    input:
        rules.align__recalibrate.input,
