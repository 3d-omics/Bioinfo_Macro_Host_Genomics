include: "reference/recompress.smk"


rule reference__all:
    """Re-bgzip the reference genome and known variants"""
    input:
        rules.reference__recompress__all.input,
