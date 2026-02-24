📁 Recommended directory structure

Keep your files organized:

nkhan3/
│
├── 1st_LongRead_Nanopore_Run/   # raw FASTQ files
│   reads.fastq
│
├── assembly/                     # PBS job scripts
│   flye_job.pbs
│
└── notebooks/
    nanopore_pipeline.ipynb

🟢 Step 1 — Set paths in Jupyter

# Set paths for raw reads, output, and job script

raw_reads = "/rds/general/user/nkhan3/home/1st_LongRead_Nanopore_Run/reads.fastq"

assembly_dir = "/rds/general/user/nkhan3/home/assembly/flye_output"

pbs_script = "/rds/general/user/nkhan3/home/assembly/flye_job.pbs"

This makes the notebook flexible — change the paths once, and all cells will work.

🟢 Step 2 — Submit PBS job from notebook

# Submit job

import subprocess

submit = subprocess.run(["qsub", pbs_script], capture_output=True, text=True)

job_id = submit.stdout.strip()

print("Submitted job ID:", job_id)

✅ This submits your Flye job exactly like from terminal.

✅ You’ll see the PBS job ID returned.

🟢 Step 3 — Monitor PBS job

# Function to check job status


def check_job_status(user="nkhan3"):

  import subprocess
  
  result = subprocess.run(["qstat", "-u", user], capture_output=True, text=True)
    
  print(result.stdout)
    

# Example usage

check_job_status()

You can run this cell repeatedly to see if the job is still running.

🟢 Step 4 — Wait for completion (optional)

You can make it semi-automatic:

import time

import subprocess

job_finished = False

while not job_finished:

  result = subprocess.run(["qstat", "-u", "nkhan3"], capture_output=True, text=True)
  
  if job_id not in result.stdout:
      
  job_finished = True

  print("Job finished!")
   
  else:
       
  print("Job still running...")
       
  time.sleep(300)  # check every 5 minutes

This way, your notebook can “pause” until the assembly finishes.

🟢 Step 5 — Inspect assembly output

Once the job is done:

import os

print("Files in assembly output folder:")

os.listdir(assembly_dir)

You should see:

assembly.fasta

assembly_info.txt

flye.log

🟢 Step 6 — Parse assembly FASTA

Example: count number of contigs and lengths

from Bio import SeqIO

fasta_file = os.path.join(assembly_dir, "assembly.fasta")

contig_lengths = []

for record in SeqIO.parse(fasta_file, "fasta"):

  contig_lengths.append(len(record.seq))

print("Number of contigs:", len(contig_lengths))

print("Contig lengths (bp):", contig_lengths[:10])  # first 10 contigs

Uses Biopython — you can install in S4 if not already:

conda activate S4

conda install biopython -y

🟢 Step 7 — Plot contig size distribution

import matplotlib.pyplot as plt

plt.figure(figsize=(10,6))

plt.hist(contig_lengths, bins=30, color='skyblue')

plt.title("Flye Assembly Contig Size Distribution")

plt.xlabel("Contig length (bp)")

plt.ylabel("Count")

plt.show()


✅ Summary of what this notebook does

Submits your PBS Flye job from notebook

Monitors job progress

Checks assembly output

Parses assembly.fasta

Generates assembly QC plots

All heavy computation runs on the cluster, login node stays safe, and your notebook serves as pipeline + documentation + visualization.
