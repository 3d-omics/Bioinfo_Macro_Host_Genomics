include: "swaps/somalier.smk"


rule swaps__all:
    input:
        rules.swaps__somalier__all.input,
