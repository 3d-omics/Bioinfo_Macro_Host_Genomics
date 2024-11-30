include: "bwamem2_functions.smk"


rule align__bwamem2__index:
    """Build genome index with bwa"""
    input:
        reference=REFERENCE / f"{HOST_NAME}.fa.gz",
    output:
        multiext(
            str(INDEX / HOST_NAME), ".amb", ".bwt.2bit.64", ".pac", ".0123", ".ann"
        ),
    log:
        INDEX / f"{HOST_NAME}.log",
    cache: "omit-software"
    threads: 8
    resources:
        mem_mb=64 * 1024,
        runtime=24 * 60,
    wrapper:
        "v5.2.1/bio/bwa-mem2/index"


rule align__bwamem2__index__all:
    input:
        [
            INDEX / f"{HOST_NAME}.{extension}"
            for extension in ["amb", "bwt.2bit.64", "pac", "0123", "ann"]
        ],


rule align__bwamem2__map:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        reads=[
            READS / "{sample_id}.{library_id}_1.fq.gz",
            READS / "{sample_id}.{library_id}_2.fq.gz",
        ],
        idx=multiext(
            str(INDEX / f"{HOST_NAME}"),
            ".amb",
            ".bwt.2bit.64",
            ".pac",
            ".0123",
            ".ann",
        ),
    output:
        MAP / "{sample_id}.{library_id}.cram",
    log:
        MAP / "{sample_id}.{library_id}.log",
    params:
        extra=lambda w: f"-R '{compose_read_group_header(w)}'",
        sort="samtools",
        sort_order="coordinate",
        sort_extra="",
    group:
        "align_{sample_id}"
    threads: 24
    resources:
        mem_mb=64 * 1024,
        runtime=24 * 60,
    wrapper:
        "v5.2.1/bio/bwa-mem2/mem"


rule align__bwamem2__map__all:
    """Collect the results of `bowtie2_map_one` for all libraries"""
    input:
        [
            MAP / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],


rule align__bwamem2__all:
    input:
        rules.align__bwamem2__index__all.input,
        rules.align__bwamem2__map__all.input,
