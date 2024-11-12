include: "swaps/somalier.smk"
include: "swaps/multiqc.smk"


rule swaps__all:
    input:
        rules.swaps__somalier__all.input,
        rules.swaps__multiqc__all.input,
