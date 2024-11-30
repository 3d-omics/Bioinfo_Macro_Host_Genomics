rule swaps__multiqc:
    input:
        SOMALIER / "relate.pairs.tsv",
        SOMALIER / "relate.samples.tsv",
    output:
        html=RESULTS / "swaps.html",
        zip=RESULTS / "swaps.zip",
    log:
        RESULTS / "swaps.log",
    params:
        extra="--title swaps --force",
    wrapper:
        "v5.1.0/bio/multiqc"


rule swaps__multiqc__all:
    """Collect all per step reports for the pipeline"""
    input:
        rules.swaps__multiqc.output,
