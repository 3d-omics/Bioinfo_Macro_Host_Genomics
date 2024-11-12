rule swaps__multiqc:
    input:
        SOMALIER / "relate.pairs.tsv",
        SOMALIER / "relate.samples.tsv",
    output:
        html=RESULTS / "swaps.html",
    log:
        RESULTS / "swaps.log",
    conda:
        "../../environments/multiqc.yml"
    params:
        dir=RESULTS,
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
        rules.report__step__swaps.output,
