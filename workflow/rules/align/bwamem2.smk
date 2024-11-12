include: "bwamem2_functions.smk"


rule align__bwamem2__index:
    """Build genome index with bwa"""
    input:
        reference=REFERENCE / "genome.fa.gz",
    output:
        multiext(f"{INDEX}/genome.", "amb", "bwt.2bit.64", "pac", "0123", "ann"),
    log:
        INDEX / "build.log",
    conda:
        "../../environments/bwamem2.yml"
    params:
        output_path=str(INDEX / "genome"),
        extra=params["align"]["bwamem2"]["extra"],
    cache: "omit-software"
    shell:
        """
        bwa-mem2 index \
            -p {params.output_path} \
            {input.reference} \
            {params.extra} \
        2> {log} 1>&2
        """


rule align__bwamem2__index__all:
    input:
        rules.align__bwamem2__index.output,


rule align__bwamem2__map:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=READS / "{sample_id}.{library_id}_2.fq.gz",
        idx=multiext(f"{INDEX}/genome", ".amb", ".bwt.2bit.64", ".pac", ".0123", ".ann"),
        reference=REFERENCE / "genome.fa.gz",
        fai=REFERENCE / "genome.fa.gz.fai",
        gzi=REFERENCE / "genome.fa.gz.gzi",
    output:
        cram=MAP / "{sample_id}.{library_id}.cram",
    log:
        MAP / "{sample_id}.{library_id}.log",
    conda:
        "../../environments/bwamem2.yml"
    params:
        index_prefix=str(INDEX / "genome"),
        extra=params["align"]["bwamem2"]["extra"],
        read_group_header=compose_read_group_header,
    group:
        "align_{sample_id}"
    shell:
        """
        (bwa-mem2 mem \
            -t {threads} \
            -R '{params.read_group_header}' \
            {params.index_prefix} \
            {params.extra} \
            {input.forward_} \
            {input.reverse_} \
        | samtools sort \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


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
