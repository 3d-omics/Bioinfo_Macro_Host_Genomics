def get_crams_for_mark_duplicates(wildcards):
    """Get all the bams for mark duplicates. They will be joined by sample"""
    libraries = samples[
        samples["sample_id"] == wildcards.sample_id
    ].library_id.values.tolist()
    crams = [
        MAP / f"{wildcards.sample_id}.{library_id}.cram" for library_id in libraries
    ]
    return crams


def compose_input_line_for_mark_duplicates(wildcards):
    """Compose input line for mark duplicates"""
    crams = get_crams_for_mark_duplicates(wildcards)
    return " ".join(f"--INPUT {cram}" for cram in crams)
