def compose_merge_vcfs_input_line(wildcards):
    files = [CALL / f"{region}.vcf.gz" for region in REGIONS]
    string = " ".join(f"--INPUT {file}" for file in files)
    return string
