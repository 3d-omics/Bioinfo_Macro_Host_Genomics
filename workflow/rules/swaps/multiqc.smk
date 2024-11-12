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


rule swaps__multiqc__all:
    """Collect all per step reports for the pipeline"""
    input:
        rules.swaps__multiqc.output,
