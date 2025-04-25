# sample1.Snakefile: FASTQ Data Download and Counting Example
# ----------------------------------------------------------
# This Snakefile shows how to download a FASTQ file, unzip it,
# and count the number of reads (headers starting with '@').

rule all:
    input:
        "sample/SRR000001.fastq.gz",
        "sample/SRR000001.fastq",
        "results/SRR000001.seq.count"

# 1) Download compressed FASTQ
rule download_reads:
    output:
        "sample/SRR000001.fastq.gz"
    shell:
        "wget -O {output} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/000/SRR000001/SRR000001.fastq.gz"

# 2) Unzip to FASTQ
rule unzip_reads:
    input:
        "sample/SRR000001.fastq.gz"
    output:
        "sample/SRR000001.fastq"
    shell:
        "gunzip -c {input} > {output}"

# 3) Count reads in FASTQ
rule count_sequences:
    input:
        "sample/SRR000001.fastq"
    output:
        "results/SRR000001.seq.count"
    shell:
        "grep -c '^@' {input} > {output}" 