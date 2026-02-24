# 🧠 Part 1 — What IGV actually needs

IGV does NOT read raw sequencing files directly.

Input you might have vs What IGV actually needs

FASTQ (raw reads)❌ Not directly usable

FASTA (contigs)	❌ Not for read viewing

SAM	⚠ Too big, slow

BAM	✅ Yes

BAM + BAI index	✅ REQUIRED

bigWig	✅ For coverage tracks

# So the standard pipeline is:

pipeline = FASTQ → Align → SAM → BAM → sorted BAM → index BAM → IGV

If input = .ab1 → convert → FASTA/FASTQ → align

If input = .seq → convert → FASTA → align

If input = .fas → align directly

If input = .fastq → align directly


# file conversions 
convert AB1 to FASTQ (better if you want qualities):

    def ab1_to_fastq(ab1_file, fastq_out):

    record = SeqIO.read(ab1_file, "abi")
    
    SeqIO.write(record, fastq_out, "fastq")
    
    print("Saved:", fastq_out)

# 🏗 Part 2 — Industry standard tools

In real labs:

Alignment: bwa or hisat2 or bowtie2

Format handling: samtools

Coverage: deeptools or bedtools

We will use:

bwa

samtools

Because:

They are universal

They are IGV-compatible

They are production-grade

# 🛠 Part 3 — First: Make sure tools are installed

In your bioinformatics conda environment:

conda install bwa samtools

Test:

bwa
samtools

# Continue setup 

⭐Step 6 — Install the Jupyter kernel (so VS Code can see it)

    python -m ipykernel install --user --name bioinfo_env

⭐ Step 7 — Open VS Code through the same terminal

    code .

Now VS Code opens in the correct environment.

⭐ Step 8 — Create a Jupyter Notebook

In VS Code:

File → New File

Save as: analysis.ipynb

⭐ Step 9 — Select your kernel

Top right → Select Kernel → choose:

bioinfo_env

⭐ Step 10 — Test example

      import pandas as pd  
      
      import seaborn as sns 
      
    import matplotlib.pyplot as plt 
    
      print("Jupyter working in Miniconda + conda-forge!")

If it prints → you’re fully set up.

⭐ Step 11 - Coming Back to the files
when re=opening the notebook, select the correct kernel.

Python will offload packages - create a .txt file called requirements.txt

List in the txt file, the packages you need.

At the very beginneing of your notebook, open a cell, and run the following code to install the packages your script depeends on

    %pip install -r requirements.txt]


# 📓 Part 4 — Jupyter Notebook Reusable Pipeline Script

You can put this directly in a Jupyter notebook cell.

✅ Cell 1 — Python helpers to run shell commands: 

          import subprocess
          
        import pathlib
        
    def run(cmd):
    
    print("Running:", cmd)
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
    
        print("ERROR:")
        
        print(result.stderr)
        
        raise RuntimeError("Command failed")
        
    else:
    
        print(result.stdout)

🧠 What this does:

This lets Python run terminal commands like bwa and samtools.

✅ Cell 2 — The reusable IGV preparation function

create a function to prepare sequencing data for visualization in IGV. 
This function aligns FASTQ files to a reference genome and generates sorted BAM files 
reference_fasta: path to the reference genome in FASTA format
fastq1: path to the first FASTQ file (or single-end FASTQ file)
fastq2: path to the second FASTQ file (for paired-end data), default is None
out_prefix: prefix for output files, default is "sample"
out_dir: directory to store output files, default is "results"
pathlib.Path is used to handle file paths and output directories
run is used to execute shell commands
each step is executed using command-line tools like BWA and Samtools

    def prepare_for_igv(reference_fasta, fastq1, fastq2=None, out_prefix="sample", out_dir="results"):
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    reference_fasta = pathlib.Path(reference_fasta)
    run(f"bwa index {reference_fasta}")

# Align reads to the reference genome

Determine if single-end or paired-end and run BWA accordingly#
single-end means only one FASTQ file is provided
paired-end means two FASTQ files are provided
if more than 2 FASTQ files are provided, raise an error
    
    sam_file = out_dir / f"{out_prefix}.sam"

    if fastq2 is None:
        # single-end
        run(f"bwa mem {reference_fasta} {fastq1} > {sam_file}")
    else:
        # paired-end
        run(f"bwa mem {reference_fasta} {fastq1} {fastq2} > {sam_file}")


# Convert SAM to BAM
BAM is a binary format that is more efficient for storage and processing than SAM

        bam_file = out_dir / f"{out_prefix}.bam"
        run(f"samtools view -bS {sam_file} > {bam_file}")
    
# Sort and index the BAM file for IGV visualization. Sorting means arranging the alignments in order of their genomic coordinates
    
    sorted_bam = out_dir / f"{out_prefix}.sorted.bam"

    # Sort and index the BAM file for IGV visualization
    run(f"samtools sort {bam_file} -o {sorted_bam}")
    run(f"samtools index {sorted_bam}")]
    
 Indexing creates a companion file that allows IGV to quickly access specific regions of the BAM file
print confirmation messages with the paths to the generated files

    print("✅ IGV-ready files created:")
    print("Sorted BAM file and index:")
    
sorted_bam means the path to the sorted BAM file

    print(sorted_bam)
    
 bam.bai is the index file for the sorted BAM filw
 
      print(sorted_bam.with_suffix(".bam.bai"))
      
 sorted_bam is returned for further use if needed. it meeans the function caller can access the path to the sorted BAM file
 
    return sorted_bam
   
What this function does (in plain English)

Makes sure output folder exists

Indexes reference genome

Aligns reads with BWA

Converts SAM → BAM

Sorts BAM

Indexes BAM

Produces:

sample.sorted.bam
sample.sorted.bam.bai

These two files are what IGV needs.

# ▶️ Part 5 — How to use it

In another notebook cell:

    prepare_for_igv(
    reference_fasta="reference/ecoli.fa",
    fastq1="raw_data/SRR2584863_1.fastq",
    fastq2="raw_data/SRR2584863_2.fastq",
    out_prefix="ecoli_test",
    out_dir="results")

# 🧪 Part 6 — How to verify success

In terminal or notebook:

    ls results

You should see:

    ecoli_test.sorted.bam
    ecoli_test.sorted.bam.bai`

# 🧬 Part 7 — How to load into IGV

Open IGV

Load reference genome (same FASTA)

Load:

    ecoli_test.sorted.bam


IGV automatically uses the .bai index
