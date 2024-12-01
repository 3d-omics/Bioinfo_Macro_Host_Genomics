include: "annotate/vep.smk"
include: "annotate/multiqc.smk"


rule annotate__all:
    input:
        rules.annotate__vep__all.input,
        rules.annotate__multiqc__all.input,
