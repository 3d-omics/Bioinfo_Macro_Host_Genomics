rule helpers__samtools__stats_cram:
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
        reference=REFERENCE / "genome.fa.gz",
        fai=REFERENCE / "genome.fa.gz.fai",
        gzi=REFERENCE / "genome.fa.gz.gzi",
    output:
        "{prefix}.stats",
    log:
        "{prefix}.stats",
    conda:
        "../../environments/samtools.yml"
    shell:
        """
        samtools stats \
            --reference {input.reference} \
            {input.cram} \
        > {output} \
        2> {log}
        """
