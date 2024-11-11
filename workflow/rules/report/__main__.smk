include: "__functions__.smk"
include: "step.smk"


rule report:
    input:
        rules.report__step.input,
