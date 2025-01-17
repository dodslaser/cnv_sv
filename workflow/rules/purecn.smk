__author__ = "Erik Demitz-Helin"
__copyright__ = "Copyright 2022, Erik Demitz-Helin"
__email__ = "erik.demitz-helin@gu.se"
__license__ = "GPL-3"


rule purecn_coverage:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
    output:
        temp(
            expand(
                "cnv_sv/purecn_coverage/{{sample}}{ext}",
                ext=[
                    "_coverage.txt.gz"
                    "_coverage.loess.txt.gz"
                    "_coverage.loess.png"
                    "_coverage.loess_qc.txt"
                ]
            )
        ),
    params:
        intervals=config.get("purecn_coverage", {}).get("intervals", ""),
        extra=config.get("purecn_coverage", {}).get("extra", ""),
    log:
        "cnv_sv/purecn_coverage/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/purecn_coverage/{sample}_{type}.output.benchmark.tsv",
            config.get("purecn_coverage", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn_coverage", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("purecn_coverage", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn_coverage", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("purecn_coverage", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn_coverage", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn_coverage", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("purecn_coverage", {}).get("container", config["default_container"])
    conda:
        "../envs/purecn_coverage.yaml"
    message:
        "{rule}: Calculate coverage for {wildcards.sample}"
    shell:
        "(Rscript $PURECN/Coverage.R "
        "--out-dir=cnv_sv/purecn_coverage "
        "--bam={input.bam} "
        "--intervals={params.intervals} "
        "--cores={threads} "
        "{params.extra}) &> {log}"


rule purecn:
    input:
        unpack(get_purecn_inputs),
        vcf="snv_indels/gatk_mutect2/{sample}_T.merged.softfiltered.vcf.gz",
        tbi="snv_indels/gatk_mutect2/{sample}_T.merged.softfiltered.vcf.gz.tbi",
    output:
        temp(
            expand(
                "cnv_sv/purecn/{{sample}}{ext}",
                ext=[
                    ".rds"
                    ".pdf"
                    "_dnacopy.seg"
                    "_chromosomes.pdf"
                    "_segmentation.pdf"
                    "_local_optima.pdf"
                    "_variants.csv"
                    "_loh.csv"
                ]
            )
        ),
    params:
        fun_segmentation=config.get("purecn", {}).get("fun_segmentation", ""),
        genome=config.get("purecn", {}).get("genome", ""),
        interval_padding=config.get("purecn", {}).get("interval_padding", ""),
        extra=get_purecn_extra,
    log:
        "cnv_sv/purecn/{sample}_T.output.log",
    benchmark:
        repeat(
            "cnv_sv/purecn/{sample}_T.output.benchmark.tsv",
            config.get("purecn", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("purecn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("purecn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("purecn", {}).get("container", config["default_container"])
    conda:
        "../envs/purecn.yaml"
    message:
        "{rule}: Quantify purity/ploidy for {wildcards.sample}"
    shell:
        "(Rscript $PURECN/PureCN.R "
        "--tumor={input.tumor} "
        "--vcf={input.vcf} "
        "--genome={params.genome} "
        "--fun-segmentation={params.fun_segmentation} "
        "--interval-padding={params.interval_padding} "
        "--sampleid={wildcards.sample} "
        "--out=cnv_sv/purecn/{wildcards.sample}_T "
        "--force --seed=1337 "
        "{params.extra}) &> {log}"
