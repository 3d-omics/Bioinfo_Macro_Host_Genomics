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


rule annotate__vep__downlaod_plugins:
    output:
        directory(VEP / "plugins"),
    params:
        release=100,
    wrapper:
        "v5.2.1/bio/vep/plugins"


rule annotate__vep:
    input:
        calls=VEP / "{sample}.vcf",
        fasta=REFERENCE / f"{HOST_NAME}.fa.gz",
        gff=REFERENCE / f"{HOST_NAME}.gtf.gz",
        gtf_tbi=REFERENCE / f"{HOST_NAME}.gtf.gz.tbi",
        plugins=VEP / "plugins",
    output:
        calls=VEP / "{sample}.vcf.gz",
        stats=VEP / "{sample}.vep.html",
    log:
        VEP / "{sample}.log",
    params:
        extra="--buffer_size 500",
        plugins=[],
    resources:
        mem_mb=16 * 1024,
        runtime=8 * 60,
    wrapper:
        "v5.2.1/bio/vep/annotate"


rule annotate__vep__all:
    input:
        [VEP / f"{sample}.vcf.gz" for sample in SAMPLES],
