rule sequali_qc:
    input:
        sample=sequali_qc_reads,
    output:
        html=f"<results>/qc/{{qc_stage}}/{{group}}_{{sample,{ALL_ROW_SAMPLE_PATTERN}}}.sequali.html",
        json="<results>/qc/{qc_stage}/{group}_{sample}.sequali.json",
    log:
        "<logs>/sequali/{qc_stage}/{group}_{sample}.sequali.log",
    wildcard_constraints:
        qc_stage="raw|processed",
    threads: 2
    params:
        extra=SEQUALI_EXTRA,
    message:
        """--- Running Sequali QC on {wildcards.qc_stage} reads."""
    wrapper:
        "v9.13.0/bio/sequali"
