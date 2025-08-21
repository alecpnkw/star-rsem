rule rsem_prepare_reference:
    input:
        fasta = config["genome"]["fasta"],
        gtf = config["genome"]["annotation"]
    output:
        "results/rsem-index/{genome}/{genome}.idx.fa"
    params:
        prefix = lambda wc, output: dirname(output[0]) + "/{0}".format(wc.genome)
    threads:
        8
    #conda: 
    #    "../envs/rsem.yaml"
    envmodules:
        "rsem/1.3.0"
    resources:
        walltime = 360
        mem_mb = 16000
    shell:
        """
        rsem-prepare-reference \
        --gtf {input.gtf} \
        --num-threads {threads} \
        {input.fasta} \
        {params.prefix}
        """

rule rsem_calculate_expression:
    input:
        genome = "results/rsem-index/{genome}/{genome}.idx.fa",
        transcripts = "results/star-pe/{sample}_{genome}/Aligned.toTranscriptome.out.bam"
    output:
        multiext("results/rsem/{sample}_{genome}/{sample}_{genome}", ".isoforms.results", ".genes.results", ".transcript.bam")
    params:
        output_prefix = lambda wc, output: dirname(output[0]) + "/{0}_{1}".format(wc.sample, wc.genome),
        genome_prefix = lambda wc, input: dirname(input["genome"]) + "/" + ".".join(basename(input["genome"]).split('.')[:-2]),
        threads = 16
    threads:
        8
    resources:
        walltime = 720
        mem_mb = 16000 
    #conda:
    #    "envs/rsem.yaml"
    envmodules:
        "rsem/1.3.0"
    shell:
        """
        rsem-calculate-expression \
        --alignments \
        --strandedness reverse \
        --paired-end \
        --num-threads {params.threads} \
        {input.transcripts} \
        {params.genome_prefix} \
        {params.output_prefix}
        """
        # Check strandedness of library, perhaps push to config!

rule gather_rsem_transcripts:
    input:
        expand("results/rsem/{sample}_{{genome}}/{sample}_{{genome}}.isoforms.results", 
        sample=samples.index
        )
    output:
        "results/rsem-combined/rsem_transcript_{genome}.tsv"
    params:
        samples = samples.index,
        value = "expected_count",
        key = "transcript_id"
    resources:
        walltime = 30
        mem_mb = 5000
    script:
        "../scripts/gather-rsem.py"

rule gather_rsem_genes:
    input:
        expand("results/rsem/{sample}_{{genome}}/{sample}_{{genome}}.genes.results", 
        sample=samples.index
        )
    output:
        "results/rsem-combined/rsem_gene_{genome}.tsv"
    params:
        samples = samples.index,
        value = "expected_count",
        key = "gene_id"
    resources:
        walltime = 30
        mem_mb = 5000
    script:
        "../scripts/gather-rsem.py"
                