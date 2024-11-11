# include: "annotate/snpeff.smk"
include: "annotate/vep.smk"
include: "annotate/multiqc.smk"


rule annotate__all:
    input:
        # rules.annotate__snpeff.input,
        rules.annotate__vep__all.input,
        rules.annotate__multiqc__all.input,
