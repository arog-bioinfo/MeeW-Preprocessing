# ==================================
# ===   short SE preprocessing   ===
# ==================================

if HAS_SHORT_SE:

    rule process_short_se_fastp:
        input:
            unpack(short_se_fastp_inputs),
        output:
            trimmed=f"<results>/reads/{{group}}_{{sample,{SHORT_SE_ROW_SAMPLE_PATTERN}}}.fastq",
            html="<results>/qc/processed/{group}_{sample}.fastp.html",
            json="<results>/qc/processed/{group}_{sample}.fastp.json",
        log:
            "<logs>/fastp/{group}_{sample}.fastp.log",
        threads: 2
        params:
            adapters=FASTP_ADAPTERS,
            extra=compose_fastp_extra,
        message:
            """--- Processing short-read SE data with fastp."""
        wrapper:
            "v9.13.0/bio/fastp"


# ==================================
# ===   short PE preprocessing   ===
# ==================================

if HAS_SHORT_PE:

    rule process_short_pe_fastp:
        input:
            unpack(short_pe_fastp_inputs),
        output:
            trimmed=[
                f"<results>/reads/{{group}}_{{sample,{SHORT_PE_ROW_SAMPLE_PATTERN}}}_R1.fastq",
                "<results>/reads/{group}_{sample}_R2.fastq",
            ],
            html="<results>/qc/processed/{group}_{sample}.fastp.html",
            json="<results>/qc/processed/{group}_{sample}.fastp.json",
        log:
            "<logs>/fastp/{group}_{sample}.fastp.log",
        threads: 2
        params:
            adapters=FASTP_ADAPTERS,
            extra=compose_fastp_extra,
        message:
            """--- Processing short-read PE data with fastp."""
        wrapper:
            "v9.13.0/bio/fastp"


# ==================================
# ===     long preprocessing     ===
# ==================================

if HAS_LONG:

    rule process_long_porechop_abi:
        input:
            unpack(long_porechop_inputs),
        output:
            fq=temp(
                f"<results>/reads/{{group}}_{{sample,{LONG_ROW_SAMPLE_PATTERN}}}.porechop.fastq"
            ),
        log:
            "<logs>/porechop_abi/{group}_{sample}.porechop_abi.log",
        threads: 2
        params:
            extra=PORECHOP_ABI_EXTRA,
            ab_initio=False,
        message:
            """--- Processing long-read data with porechop_abi."""
        wrapper:
            "v9.13.0/bio/porechop_abi"

    rule process_long_chopper:
        input:
            unpack(chopper_inputs),
        output:
            out_reads=f"<results>/reads/{{group}}_{{sample,{LONG_ROW_SAMPLE_PATTERN}}}.fastq",
        log:
            "<logs>/chopper/{group}_{sample}.chopper.log",
        threads: 2
        params:
            extra=CHOPPER_EXTRA,
        message:
            """--- Processing long-read data with chopper."""
        wrapper:
            "v9.13.0/bio/chopper"
