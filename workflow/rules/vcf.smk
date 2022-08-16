rule vcf_fix_bq:
    input:
        vcf="{file}.vcf.gz",
    output:
        vcf=temp("{file}.fix_bq.vcf"),
    benchmark:
        repeat(
            "{file}.fix_bq.benchmark.tsv",
            config.get("vcf_fix_bq", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("vcf_fix_bq", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("vcf_fix_bq", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vcf_fix_bq", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("vcf_fix_bq", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vcf_fix_bq", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vcf_fix_bq", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("vcf_fix_bq", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: Convert Mutect2 MBQ to QUAL in {output}"
    shell:
        "(gzip -dkc {input.vcf} | sed -E 's/MBQ=[0-9]+,([0-9]+)/QUAL=\\1/') > {output}"