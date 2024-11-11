rule annotate__multiqc:
    """Collect all reports for the snpeff step"""
    input:
        [VEP / f"{sample}.vep.html" for sample in SAMPLES],
    output:
        html=STEP / "annotate.html",
    log:
        STEP / "annotate.log",
    conda:
        "../../environments/multiqc.yml"
    params:
        dir=STEP,
    shell:
        """
        multiqc \
            --title annotate \
            --force \
            --filename annotate \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule annotate__multiqc__all:
    input:
        rules.annotate__multiqc.output,
