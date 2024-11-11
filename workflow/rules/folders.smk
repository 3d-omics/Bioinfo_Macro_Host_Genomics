RESULTS = Path("results/")

REFERENCE = Path("results/reference/")

ALIGN = RESULTS / "align"
READS = ALIGN / "reads"
INDEX = ALIGN / "index"
MAP = ALIGN / "map"
MARK_DUPLICATES = ALIGN / "mark_duplicates"
BCFTOOLS = ALIGN / "bcftools"
RECALIBRATE = ALIGN / "recalibrate"


VARIANTS = RESULTS / "variants"
CALL = VARIANTS / "call"
GENOTYPE = VARIANTS / "genotype"
FILTER = VARIANTS / "filter"
POSTERIORS = VARIANTS / "posteriors"

ANNOTATE = RESULTS / "annotate"
SNPEFF_DB = Path("resources/snpeff/")
SNPEFF = ANNOTATE / "snpeff"
VEP = ANNOTATE / "vep"

REPORT = Path("reports/")
STEP = REPORT / "step"

SWAPS = Path("results/swaps/")
SOMALIER = SWAPS / "somalier"
