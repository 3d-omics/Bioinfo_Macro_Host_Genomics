rule annotate__vep__tmp_vcf:
    input:
        FILTER / "all.filtered.vcf.gz",
    output:
        temp(VEP / "{sample}.vcf"),
    log:
        VEP / "{sample}.tmp_vcf.log",
    params:
        sample=lambda w: f"--samples {w.sample} --trim-alt-alleles",
    wrapper:
        "v5.2.1/bio/bcftools/view"


rule annotate__vep:
    input:
        vcf=VEP / "{sample}.vcf",
        fa=REFERENCE / f"{HOST_NAME}.fa.gz",
        gtf=REFERENCE / f"{HOST_NAME}.gtf.gz",
        gtf_tbi=REFERENCE / f"{HOST_NAME}.gtf.gz.tbi",
    output:
        vcf=VEP / "{sample}.vcf.gz",
        html=VEP / "{sample}.vep.html",
    log:
        VEP / "{sample}.log",
    conda:
        "../../environments/vep.yml"
    params:
        sample=lambda w: w.sample,
    shell:
        """
        vep \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --fasta {input.fa} \
            --gtf {input.gtf} \
            --compress_output bgzip \
            --stats_file {output.html} \
            --buffer_size 500 \
        2>> {log} 1>&2
        """


rule annotate__vep__all:
    input:
        [VEP / f"{sample}.vcf.gz" for sample in SAMPLES],
