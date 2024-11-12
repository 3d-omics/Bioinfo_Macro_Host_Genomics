rule align__multiqc:
    """Collect all reports for the align submodule step"""
    input:
        reads=[
            READS / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
        bwamem2=[
            MAP / f"{sample_id}.{library_id}.stats"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        markduplicates=[MARK_DUPLICATES / f"{sample_id}.stats" for sample_id in SAMPLES],
        recalibrate=[
            RECALIBRATE / f"{sample_id}.stats"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        html=RESULTS / "align.html",
    log:
        RESULTS / "align.log",
    conda:
        "../../environments/multiqc.yml"
    params:
        dir=RESULTS,
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
        RESULTS / "align.html",
