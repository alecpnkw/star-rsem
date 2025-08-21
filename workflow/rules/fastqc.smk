def input_fastq(wc):
    if wc.read == "R1":
        return samples.loc[wc.sample, "R1_fastq"]
    elif wc.read == "R2":
        return samples.loc[wc.sample, "R2_fastq"]

# using local script over wrapper for lsf compatibility 
rule fastqc:
    input:
        input_fastq
    output:
        html="results/fastqc/{sample}_{read}.html",
        zip="results/fastqc/{sample}_{read}_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{sample}_{read}.log"
    conda:
        "../envs/fastqc.yaml"
    envmodules:
        "fastqc/0.11.8"
    threads: 1
    resources:
        walltime = 60
        mem_mb = 4000
    script:
        "../scripts/fastqc-wrapper.py"
