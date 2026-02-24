1️⃣ Variant Calling Pipeline (for DNA sequencing)

If your goal is to find differences from a reference genome (like SNPs or indels), here’s the next workflow:

Step 1: Sort & Index BAM (if not done yet)

samtools sort results/ecoli.bam -o results/ecoli.sorted.bam
samtools index results/ecoli.sorted.bam


samtools sort → sort reads by genomic coordinates.

samtools index → creates an index file so tools can quickly access regions.

Step 2: Variant Calling with bcftools

bcftools mpileup -f reference/ecoli.fa results/ecoli.sorted.bam | \
bcftools call -mv -Oz -o results/ecoli_variants.vcf.gz


Explanation:

mpileup → looks at coverage at each position.

call -mv → calls SNPs (-v) and indels (-m).

-Oz → outputs compressed VCF (.vcf.gz).

Step 3: Index VCF

bcftools index results/ecoli_variants.vcf.gz

Step 4: Explore Variants

bcftools view results/ecoli_variants.vcf.gz | head


Shows first few variants in terminal.

For full analysis, you can import the VCF into Python or R for plotting.

2️⃣ Gene Expression / RNA-Seq Pipeline (for RNA sequencing)

If the dataset were RNA-Seq, you’d usually:

Align reads (already done with HISAT2/STAR or BWA).

Count reads per gene (e.g., using htseq-count or featureCounts).

Example with HTSeq:

htseq-count -f bam -r pos -s no results/ecoli.sorted.bam reference/ecoli.gff > results/gene_counts.txt


Explanation:

-f bam → input format is BAM.

-r pos → sorted by position.

-s no → unstranded library.

Then you can load gene_counts.txt into Python or R for differential expression analysis with DESeq2 or edgeR.

3️⃣ Optional: Summarize Results in a Notebook

Open VS Code or Jupyter notebook.

Load results:

import pandas as pd

# Example for RNA-seq counts
counts = pd.read_csv("results/gene_counts.txt", sep="\t", header=None)
counts.head()


For DNA variants, use pysam or vcf Python libraries to parse .vcf.gz and plot variant distributions.

✅ Summary

You now have a full practice pipeline:

Download test data → fasterq-dump SRR2584863

Quality control → fastqc + multiqc

Alignment → bwa mem → sorted BAM

Variant calling → bcftools mpileup + call
or RNA-seq counting → htseq-count

Explore/analyze results → Python/Jupyter or R
