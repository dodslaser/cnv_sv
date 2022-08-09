
__author__ = "Erik Demitz-Helin"
__copyright__ = "Copyright 2022, Erik Demitz-Helin"
__email__ = "erik.demitz-helin@gu.se"
__license__ = "GPL-3"



rule purecn:
    input:
        seg="cnv_sv/cnvkit_export_seg/{sample}_{type}.seg",
        regions="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr",
        vcf="snv_indels/gatk_mutect2/{sample}_{type}.merged.softfiltered.vcf.gz",
        tbi="snv_indels/gatk_mutect2/{sample}_{type}.merged.softfiltered.vcf.gz.tbi",
    output:
        temp("cnv_sv/purecn/{sample}_{type}.csv"),
    params:
        extra=config.get("purecn", {}).get("extra", ""),
    log:
        "cnv_sv/purecn/{sample}_{type}.output.log"
    benchmark:
        repeat(
            "cnv_sv/purecn/{sample}_{type}.output.benchmark.tsv", config.get("purecn", {}).get("benchmark_repeats", 1)
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
        """
        (
            Rscript $PURECN/PureCN.R \
                --sampleid={wildcards.sample}_{wildcards.type} \
                --seg-file={input.seg} \
                --tumor={input.regions} \
                --vcf={input.vcf} \
                --genome=hg19 \
                --parallel \
                --cores={threads} \
                --force \
                --seed=123 \
                --post-optimize \
                --fun-segmentation=Hclust \
                --out=cnv_sv/purecn/{wildcards.sample}_{wildcards.type} \
                {params.extra}
        ) &> {log} || (
            echo '"Sampleid","Purity","Ploidy","Sex","Contamination","Flagged","Failed","Curated","Comment"' > {output}
            echo '{wildcards.sample}_{wildcards.type},?,?,?,?,TRUE,TRUE,FALSE,PURECN_ERROR' >> {output}
        )
        """

