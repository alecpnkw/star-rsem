def input_fastq(wc):
    if wc.read == "R1":
        return samples.loc[wc.sample, "R1_fastq"]
    elif wc.read == "R2":
        return samples.loc[wc.sample, "R2_fastq"]

rule fastqc:
    input:
        input_fastq
    output:
        html="results/fastqc/{sample}_{read}.html",
        zip="results/fastqc/{sample}_{read}_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{sample}_{read}.log"
    threads: 1
    wrapper:
        "0.78.0/bio/fastqc" # use wrapper for handling of fnames
