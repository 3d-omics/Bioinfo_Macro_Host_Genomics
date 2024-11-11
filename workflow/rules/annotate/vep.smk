rule annotate__vep__tmp_vcf:
    input:
        vcf=FILTER / "all.filtered.vcf.gz",
    output:
        tmp=temp(VEP / "{sample}.tmp.vcf.gz"),
    log:
        VEP / "{sample}.tmp_vcf.log",
    conda:
        "../../environments/bcftools.yml"
    params:
        sample=lambda w: w.sample,
    shell:
        """
        bcftools view \
            --samples {params.sample} \
            --output-type z1 \
            --output-file {output.tmp} \
            --trim-alt-alleles \
            {input.vcf} \
        2> {log} 2>&1
        """


rule annotate__vep:
    input:
        vcf=VEP / "{sample}.tmp.vcf.gz",
        fa=REFERENCE / "genome.fa.gz",
        gtf=REFERENCE / "annotation.gtf.gz",
        gtf_tbi=REFERENCE / "annotation.gtf.gz.tbi",
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
