configfile: "config.yaml"
def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule bwa_map: #all rules use this format with an input and an output
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
        # "data/samples/{sample}.fastq"
        #{sample} allows you to use this rule for multiple fastq files
        #for example snakemake -np mapped_reads/{A,B}.bam the {A,B} allows you to list multiple output file names
        #touch data/samples/A.fastq updates file
    output:
        temp("mapped_reads/{sample}.bam") #in a snakemake file you can name a new directory and it will automatically know to make the new directory if the directory doesn't exist yet
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"    
    benchmark: 
        "benchmarks/{sample}.bwa.benchmark.txt" #indicates where to store benchmarking results (clock time and memory usage in MiB stored in tab delimited format) to perform multiple analyses in the same file in subsequent lines use repeat("file_name", #times)
    # threads: 8
    # wrapper:
    #   "0.15.3/bio/bwa/mem"
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"

rule samtools_sort: #sorts reads and stores a sorted version in sorted_reads directory.... order?
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")#protected makes certain that it can't be overwritten #never overwriting, new file for each step, this is a good practice, where each rule creates a new directory
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} " #this is the command that will actually be sent to the shell ie "snakemake -np sorted_reads/B.bam -O bam mapped_reads/B.bam"
        #wildcards allows you to acess the objects attribute of sample number or name
        "-O bam {input} > {output}"