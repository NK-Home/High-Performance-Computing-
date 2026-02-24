## Set up a small test dataset from NCBI SRA (Sequence Read Archive):

## 🔧 Step 1: Install the SRA Toolkit

This lets you fetch data directly from NCBI. In your bioinformatics environment:

mamba install -c bioconda sra-tools


- Test it:

fastq-dump --version

You should see the version number.

## 🔧 Step 2: Pick a Small Example Dataset

We’ll use a tiny E. coli dataset that’s often used for tutorials:

SRA Run Accession: SRR2584863
(paired-end Illumina sequencing of E. coli, just ~1MB — very small).

## 🔧 Step 3: Download FASTQ Files

Run:
fasterq-dump SRR2584863 --split-files -O raw_data/

Explanation:

fasterq-dump → improved version of fastq-dump (faster, multi-threaded).

--split-files → if the run is paired-end, you’ll get _1.fastq and _2.fastq.

-O raw_data/ → saves into your raw_data folder.

Check files:

ls -lh raw_data/

You should see something like:

SRR2584863_1.fastq

SRR2584863_2.fastq

- if it gets stuck - use this instead = prefetch first -
  
prefetch SRR2584863 -O raw_data/

prefetch downloads the SRA file locally before converting to FASTQ — often more reliable.

# Step 2: Convert to FASTQ
fasterq-dump raw_data/SRR2584863/SRR2584863.sra --split-files -O raw_data/

prefetch downloads the .sra archive.

fasterq-dump then converts it to FASTQ locally, which avoids network issues during conversion.

Check files:

ls -lh raw_data/

## 🔧 Step 4: Run Quality Control

Now you can test your pipeline:

fastqc raw_data/SRR2584863_*.fastq -o results/
multiqc results/

This reads the FASTQ files and generates .html reports and .zip summary files.

## check out the html files 

To open the HTML reports in your Windows browser:

explorer.exe results/multiqc_report.html

## 🔧 Step 5: Try Alignment

For practice, align against an E. coli reference genome.
Download the reference FASTA:

wget -qO- \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz \
  | gunzip > reference/ecoli.fa

Download the genome FASTA from NCBI (wget).

wget is a command-line tool to download files from the web.

These are options (flags) for wget:

-q = quiet mode
Suppresses extra output like download progress bars. (Keeps your terminal clean.)

-O- = write output to stdout (the terminal), not to a file
The -O option normally specifies a file name. Example:

wget -O ecoli.fa.gz https://...

But here we use a dash -, which means “send the data to standard output (stdout).”
This allows us to pipe it into another command instead of saving directly.

| gunzip

| = pipe operator.
It takes the output from the command on the left (wget) and passes it as input to the command on the right (gunzip).

gunzip = decompresses .gz files.


# Index it:

bwa index reference/ecoli.fa

BWA = Burrows-Wheeler Aligner.

It’s a popular program for aligning sequencing reads (FASTQ) to a reference genome (FASTA).

index

This tells BWA to prepare the reference genome (ecoli.fa) for fast searching.

When aligning, BWA needs to look up where small read sequences (like 100 bp Illumina reads) match inside the huge genome (millions of bases).

After bwa index, if you check your reference/ folder:

ls reference/ 

they’re required by BWA when you run bwa mem for alignment.

# Align reads:

bwa mem reference/ecoli.fa raw_data/SRR2584863_1.fastq raw_data/SRR2584863_2.fastq | samtools view -bS - > results/ecoli.bam

1. bwa mem

Runs BWA (Burrows-Wheeler Aligner) in MEM mode.

mem = Maximal Exact Matches, the algorithm optimized for Illumina-style reads (100–150 bp).

So, this is the aligner step: it tries to find where each short read belongs in the genome.

does this againt the reference genome and the SRA.files 
Each read pair comes from opposite ends of the same DNA fragment, which helps improve alignment accuracy.

| samtools view -bS -

This part takes the output from BWA and processes it with Samtools:

| = pipe → send the output of BWA to the next command instead of printing it to screen.

samtools view = converts alignment formats.

By default, BWA outputs SAM format (text-based, big file).

We convert it into BAM format (binary, compressed, faster).

Options:

-b = output BAM (not text SAM).

-S = input is SAM.

- = read input from stdin (the pipe from BWA).

5. > results/ecoli.bam

Redirects the final BAM alignment file into your results/ directory.

ecoli.bam = contains all your reads mapped to the E. coli reference genome.

# Sort & index:

samtools sort results/ecoli.bam -o results/ecoli.sorted.bam
samtools index results/ecoli.sorted.bam

## IGV (Integrative Genomics Viewer) lets you see reads aligned to the genome.

Step 1. Install IGV

If you’re on Windows (with Ubuntu/WSL):

Download IGVmfor your system.

Install and run it from Windows (not WSL).

Step 2. Copy or open files

You’ll need:

The reference genome you aligned to: reference/ecoli.fa

The BAM file: results/ecoli.sorted.bam

The BAM index: results/ecoli.sorted.bam.bai

If they are inside WSL (~/bioinformatics/...), you can access them in Windows here:

\\wsl$\Ubuntu\home\<your-username>\bioinformatics\


and then load them into IGV.

## ✅ Now you have:

QC report (multiqc_report.html).

Aligned reads (ecoli.sorted.bam).

This is the full mini workflow from raw data → QC → alignment.
