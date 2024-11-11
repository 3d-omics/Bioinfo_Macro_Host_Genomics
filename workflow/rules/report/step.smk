rule report__step__reads:
    """Collect all reports for the reads step"""
    input:
        [
            READS / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
    output:
        html=STEP / "reads.html",
    log:
        STEP / "reads.log",
    conda:
        "__environment__.yml"
    params:
        dir=STEP,
    shell:
        """
        multiqc \
            --filename reads \
            --title reads \
            --force \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__align:
    """Collect all reports for the bowtie2 step"""
    input:
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
        "__environment__.yml"
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


rule report__step__annotate:
    """Collect all reports for the snpeff step"""
    input:
        rules.annotate__vep__reports.input,
    output:
        html=STEP / "annotate.html",
    log:
        STEP / "annotate.log",
    conda:
        "__environment__.yml"
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
        "__environment__.yml"
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
        rules.report__step__reads.output,
        rules.report__step__align.output,
        rules.report__step__annotate.output,
        rules.report__step__swaps.output,
