rule annotate__multiqc:
    """Collect all reports for the VEP step"""
    input:
        [VEP / f"{sample}.vep.html" for sample in SAMPLES],
    output:
        html=RESULTS / "annotate.html",
        zip=RESULTS / "annotate.zip",
    log:
        RESULTS / "annotate.log",
    params:
        extra="--title annotate --force",
    wrapper:
        "v5.1.0/bio/multiqc"


rule annotate__multiqc__all:
    input:
        rules.annotate__multiqc.output,
