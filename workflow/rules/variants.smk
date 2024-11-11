include: "variants/call.smk"
include: "variants/genotype.smk"
include: "variants/filter.smk"


rule variants__all:
    input:
        rules.variants__filter.input,
