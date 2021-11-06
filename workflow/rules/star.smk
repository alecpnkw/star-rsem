rule star_index:
    input:
        fasta = config["genome"]["fasta"],
        gtf = config["genome"]["annotation"]
    output:
        directory("results/star-index/{genome}")
    message:
        "Building STAR index..."
    threads: 16
    #conda:
    #    "../envs/star.yaml"
    envmodules:
        "star/2.7.9a"
    shell:
        """
        STAR --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf}
        """

rule star_pe_single:
    input:
        genome = "results/star-index/{genome}",
        R1 = lambda wc: samples.loc[wc.sample, "R1_fastq"], # directs to fastq files
        R2 = lambda wc: samples.loc[wc.sample, "R2_fastq"]
    output:
        # see STAR manual for additional output files
        "results/star-pe/{sample}_{genome}/Aligned.sortedByCoord.out.bam",
        "results/star-pe/{sample}_{genome}/Aligned.toTranscriptome.out.bam",
        "results/star-pe/{sample}_{genome}/SJ.out.tab"
    threads: 12
    params:
        prefix = "results/star-pe/{sample}_{genome}/",
        threads = 24
    resources:
        mem_mb = 5000
    #conda: 
     #   "../envs/star.yaml"
    envmodules:
        "star/2.7.9a"
    shell:
        """
        STAR \
		--genomeDir {input.genome} \
		--readFilesIn {input.R1} {input.R2} \
		--readFilesCommand zcat \
		--outFileNamePrefix {params.prefix} \
		--runThreadN {params.threads} \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate
        """   