rule align__map__:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=READS / "{sample_id}.{library_id}_2.fq.gz",
        idx=multiext(f"{INDEX}/genome", ".amb", ".bwt.2bit.64", ".pac", ".0123", ".ann"),
        reference=REFERENCE / "genome.fa.gz",
    output:
        cram=MAP / "{sample_id}.{library_id}.cram",
    log:
        MAP / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        index_prefix=str(INDEX / "genome"),
        extra=params["bowtie2"]["extra"],
        samtools_mem=params["bowtie2"]["samtools"]["mem_per_thread"],
        read_group_header=compose_read_group_header,
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
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


# There is no point in separating bwa from samtools since they switch between
# aligning and sorting


rule align__map:
    """Collect the results of `bowtie2_map_one` for all libraries"""
    input:
        [
            MAP / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
