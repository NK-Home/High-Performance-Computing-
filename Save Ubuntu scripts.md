How to save your Ubuntu script
1. Open a text editor inside Ubuntu

You can use nano (simplest editor). Example:

nano my_pipeline.sh


This opens a blank file called my_pipeline.sh.

2. Add your commands

For example, paste in the commands you’ve been using:

#!/bin/bash
# A simple bioinformatics pipeline

# Step 1. Run FastQC
fastqc raw_data/SRR2584863_*.fastq -o results/

# Step 2. MultiQC summary
multiqc results/

# Step 3. Align reads with BWA
bwa mem reference/ecoli.fa raw_data/SRR2584863_1.fastq raw_data/SRR2584863_2.fastq | samtools view -bS - > results/ecoli.bam

# Step 4. Sort and index BAM
samtools sort results/ecoli.bam -o results/ecoli.sorted.bam
samtools index results/ecoli.sorted.bam

# Step 5. Get alignment stats
samtools flagstat results/ecoli.sorted.bam

3. Save the file

In nano, press:

CTRL + O → save (write out)

Enter → confirm filename

CTRL + X → exit

Now your script is saved in your current directory as my_pipeline.sh.

4. Make the script executable
chmod +x my_pipeline.sh

5. Run the script
./my_pipeline.sh


It will execute each command one after another, just like you typed them manually.

✅ Benefits

You don’t lose your work — all commands are stored.

Easy to rerun if something fails.

You can share your script with others (or future you).

Keeps your analysis reproducible.
