# sample2.Snakefile: Mapping, Sorting and Indexing Example
# -------------------------------------------------------
# This Snakefile shows how to align reads to a reference,
# sort the resulting BAM, and create an index for downstream use.

rule all:
    input:
        "sample/SRR000001.fastq",
        "alignments/SRR000001.bam",
        "alignments/SRR000001.sorted.bam",
        "alignments/SRR000001.sorted.bam.bai"

# 1) Align reads to reference (bwa)
rule map_reads:
    input:
        reads="sample/SRR000001.fastq",
        ref="data/ref.fa"
    output:
        bam="alignments/SRR000001.bam"
    shell:
        "bwa mem {input.ref} {input.reads} | samtools view -Sb - > {output.bam}"

# 2) Sort the BAM file (samtools)
rule sort_bam:
    input:
        bam="alignments/SRR000001.bam"
    output:
        sorted_bam="alignments/SRR000001.sorted.bam"
    shell:
        "samtools sort {input.bam} -o {output.sorted_bam}"

# 3) Index the sorted BAM (samtools)
rule index_bam:
    input:
        sorted_bam="alignments/SRR000001.sorted.bam"
    output:
        bai="alignments/SRR000001.sorted.bam.bai"
    shell:
        "samtools index {input.sorted_bam}" 