# Open a nano file 

nano flye_job.pbs

# Save and exit nano

Press:

CTRL + O
Enter
CTRL + X


#!/bin/bash
#PBS -N flye_assembly
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o flye_output.log

# ==========================
# Go to the directory where the job was submitted
# ==========================
cd $PBS_O_WORKDIR

echo "================================="
echo "Flye assembly job started"
date
hostname
echo "================================="

# ==========================
# Load conda and activate environment
# ==========================
source ~/miniforge3/etc/profile.d/conda.sh
conda activate S4

echo "Environment check:"
which python
which porechop
which flye
echo "================================="

# ==========================
# Step 1: Run Porechop (adapter trimming)
# ==========================
# Input: 13 GB nanopore FASTQ file
# Output: trimmed FASTQ in working directory
# Uses 16 threads
# ==========================
echo "Step 1: Running Porechop..."
porechop \
  -i ~/1st_LongRead_Nanopore_Run/reads.fastq \
  -o trimmed.fastq \
  --threads 16

# ==========================
# Step 2: Run Flye assembly
# ==========================
# Input: trimmed reads
# Output: assembly results in flye_output/
# Memory and threads optimized for 13 GB reads
# ==========================
echo "Step 2: Running Flye assembly..."
flye \
  --nano-raw trimmed.fastq \
  --out-dir flye_output \
  --genome-size 5m \
  --threads 16

# ==========================
# Step 3: Job finished
# ==========================
echo "Assembly completed!"
date
echo "Output folder: flye_output/"
echo "================================="




# Monitor progress:

qstat -u nkhan3
tail -f flye_output.log

# Monitor job

qstat -u nkhan3

# Check results

# After completion:

ls flye_output/
