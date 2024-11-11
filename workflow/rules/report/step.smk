rule report__step__annotate:
    """Collect all reports for the snpeff step"""
    input:
        rules.annotate__vep__reports.input,
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


rule report__step__swaps:
    input:
        rules.swaps__somalier__report.input,
    output:
        html=STEP / "swaps.html",
    log:
        STEP / "swaps.log",
    conda:
        "../../environments/multiqc.yml"
    params:
        dir=STEP,
    shell:
        """
        multiqc \
            --title swaps \
            --force \
            --filename swaps \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report__step__annotate.output,
        rules.report__step__swaps.output,
