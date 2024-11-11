rule align__multiqc:
    """Collect all reports for the align submodule step"""
    input:
        [
            READS / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
        [
            MAP / f"{sample_id}.{library_id}.{report}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in BAM_REPORTS
        ],
        [
            MARK_DUPLICATES / f"{sample_id}.{report}"
            for sample_id in SAMPLES
            for report in BAM_REPORTS
        ],
        [
            RECALIBRATE / f"{sample_id}.{report}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in BAM_REPORTS
        ],
    output:
        html=STEP / "align.html",
    log:
        STEP / "align.log",
    conda:
        "../../environments/multiqc.yml"
    params:
        dir=STEP,
    shell:
        """
        multiqc \
            --title align \
            --force \
            --dirs \
            --dirs-depth 1 \
            --filename align \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule align__multiqc__all:
    input:
        html=STEP / "align.html",
