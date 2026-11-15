rule infer_auto_preprocess:
    input:
        sequali_json="<results>/qc/raw/{group}_{sample}.sequali.json",
        config_schema=workflow.source_path("../schemas/config.schema.yaml"),
        sequali_schema=workflow.source_path("../schemas/sequali_short_read.schema.yaml"),
    output:
        auto_params=f"<results>/reads/{{group}}_{{sample,{SHORT_ROW_SAMPLE_PATTERN}}}.auto_params.json",
    log:
        "<logs>/auto_preprocess/{group}_{sample}.auto_params.log",
    conda:
        "../envs/auto_preprocess.yaml"
    params:
        data_type=auto_preprocess_data_type,
        auto_preprocess_config=AUTO_PREPROCESS_CONFIG,
    message:
        """--- Inferring auto preprocess parameters."""
    script:
        "../scripts/infer_preprocess_params.py"
