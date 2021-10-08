rule multiqc:
    input:
        expand("results/rsem/{sample}_{genome}/{sample}_{genome}.isoforms.results", sample=samples.index, genome = config["genome"]["build"]),
        expand("results/fastqc/{sample}_{read}.html", sample=samples.index, read = ["R1", "R2"])
    output:
        "results/multiqc/multiqc-report.html"
    shell:
        """
        multiqc \ 
        --no-data-dir \
        --filename {output}
        "results/"
        """