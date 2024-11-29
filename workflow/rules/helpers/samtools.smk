rule helpers__samtools__stats_cram:
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
        reference=REFERENCE / f"{HOST_NAME}.fa.gz",
        fai=REFERENCE / f"{HOST_NAME}.fa.gz.fai",
        gzi=REFERENCE / f"{HOST_NAME}.fa.gz.gzi",
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
