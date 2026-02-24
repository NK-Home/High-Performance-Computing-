# Introduction to Linux OS
 - linux is used because it is; free, no antivirus, avaliable programmes, doesnt slow down
 
 - Common Linux OS for bioinformatics = Cent0S and Ubuntu

# Ubuntu for Bioinformatics 
- supports tools like BLAST, STAR, Nextflow
- easy install for pacakages with APT, Conda, Docker
- Cloud integration

Bioinformatics Pipelines are run from the terminal = program to allow user to interface with shell (interperts user instructions and executes on OS)
Shell examples = BASH, CSH

# Running Linux on Windows 
- Steps to Install Ubuntu on Windows using Windows Subsystem for Linux (WSL)
  
1. Enable WSL:
Open PowerShell as an administrator (Open the Start menu, type Windows PowerShell, select Windows PowerShell, and then select Run as administrator)
- Run the command [wsl --install]. This enables WSL and install the necessary components, including a Linux kernel
- reboot the computer
  
2. Install Ubuntu:
Open the Microsoft Store, search for "Ubuntu", and install your desired version (e.g., Ubuntu 22.04 LTS). 

3. Launch Ubuntu:
After installation, you can launch Ubuntu from the Start Menu. 

4. Set up User Account:
You will be prompted to create a username and password for your Ubuntu user account.

# 🔧 Update & Secure Ubuntu

5. First, make sure your system is up to date:

sudo apt update && sudo apt upgrade -y

apt update refreshes the list of software packages.
apt upgrade installs the newest versions of your existing software.
The -y just means “yes to all prompts.”

# install some essentials:

sudo apt install build-essential curl wget git unzip htop tree -y

build-essential → compilers (needed if a tool requires building from source).

curl + wget → download data and software from the internet.

git → version control system (collaboration, reproducibility).

unzip → unpacks compressed data files (common in bioinformatics).

htop → shows running processes (important when jobs hang).

tree → shows folder structures (useful for exploring projects).

👉 Why now?
These tools form the “foundation” of any computational scientist’s toolbox. You’ll use them every day.

# 📦 Install a Package Manager (Miniforge + Mamba)

Managing bioinformatics software manually is painful (different tools need conflicting libraries). That’s why the field relies heavily on Conda (and its faster cousin, Mamba).

- Download and install Miniforge, a lightweight Conda distribution:

wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh

run yes and allow conda to make changes and run on start up 

- Restart your terminal, then:

conda install -n base -c conda-forge mamba

conda = package/environment manager.

mamba = a fast reimplementation of Conda (saves time).

if a verson already exisits then run => to update = $ conda update -n base -c conda-forge conda

👉 Why now?
Almost all bioinformatics tools are available via bioconda, which uses Conda/Mamba. Installing this first ensures every tool is just a one-line install later.

# 🧰 Create a Bioinformatics Environment & Install Core Tools

Always isolate tools into environments so they don’t conflict.

mamba create -n bioinfo -c bioconda -c conda-forge python=3.10
conda activate bioinfo

-n bioinfo → gives the environment a name.

-c bioconda -c conda-forge → tells Conda to look in the best channels for scientific software.

python=3.10 → sets a stable Python version.

- activate the enviroment - mamba activate bioinfo or mamba run -n bioinfo
  
- check if you are in the enviroment by running => conda env list
  
# install common tools:

mamba install -c conda-forge -c bioconda \
fastqc multiqc bwa samtools bedtools \
htseq bcftools seqtk hisat2 bowtie2 salmon \
trimmomatic blast mafft clustalw

-c conda-forge → provides dependencies (general science packages).

-c bioconda → provides the bioinformatics tools themselves.

Order matters: we put conda-forge before bioconda because many bioconda packages depend on conda-forge libraries.

FastQC / MultiQC → assess sequencing read quality.

BWA, HISAT2, Bowtie2 → align sequencing reads to genomes.

Samtools, Bcftools → manipulate alignments (BAM/VCF files).

Bedtools → genomic interval operations.

HTSeq → RNA-seq read counting.

Seqtk → quick sequence manipulations.

Trimmomatic → trim adapters and low-quality reads.

BLAST → search sequence databases.

MAFFT / ClustalW → multiple sequence alignments.

- After installation, check a couple of tools to confirm:
fastqc --version
samtools --version

If they output version numbers → you’re good. 🎉

👉 Why now?
These are the “daily drivers” of genomics and transcriptomics analysis. Installing them early means you can immediately start experimenting with real datasets.

##⚙️ Install Workflow Managers (for Automation)

Bioinformatics analyses often involve dozens of steps. Workflow managers make them reproducible.

Choose one (or both):

mamba install snakemake
or
mamba install nextflow


👉 Why now?
Even though you may not use them immediately, learning early will save you from repeating mistakes manually later.

## 📊 Set Up R and Python for Data Analysis

Most downstream analysis is done in R or Python.

mamba install r-base r-essentials bioconductor-deseq2 bioconductor-edger
mamba install jupyterlab pandas matplotlib seaborn scikit-learn

R + Bioconductor → specialized for genomics (DESeq2, edgeR = differential expression analysis).
Python + Jupyter → flexible for general analysis and machine learning.
Pandas / Matplotlib / Seaborn / Scikit-learn → handle data, plot results, and build predictive models.

👉 Why now?
Raw data is useless until you interpret it. These tools let you analyze and visualize results.

# 📂 Organize Your Workspace

Create a tidy directory structure:
mkdir -p ~/bioinformatics/{raw_data,results,scripts,reference,notebooks}

mkdir → makes directories (folders).

-p → tells it to create parent folders as needed (so it won’t complain if some already exist).

~ → is shorthand for your home directory.

raw_data → FASTQ/FASTA/VCF input files.

results → processed data outputs.

scripts → your Bash/Python scripts.

reference → genomes, annotation files (GTF/GFF).

notebooks → Jupyter notebooks / RMarkdown reports.

- make sure you are in the ideal directory using the code => cd ~/bioinformatics


👉 Why now?
If you don’t set up structure at the start, things quickly become chaotic. This makes your work reproducible.

# 🧪 Run Your First Mini Workflow (QC)

Let’s test the tools. Suppose you have a FASTQ file (sample.fastq.gz):

fastqc sample.fastq.gz -o results/
multiqc results/ \
&& multiqc results/ \
&& explorer.exe results/multiqc_report.html

fastqc → generates a quality report for each FASTQ file.
multiqc → combines all QC reports into one summary.
&& Operator = it takes into account the state of the previous command
explorer.exe = opens the file as its a html 

👉 Why now?
This gives you your first taste of real analysis — and lets you confirm your installation works.

# 🔄 Learn Git & Version Control

Version control ensures your scripts, notebooks, and workflows are tracked and shareable.

cd ~/bioinformatics
git init
git config --global user.name "Your Name"
git config --global user.email "you@example.com"

👉 Why now?
Science is collaborative. Git lets you work with others, track changes, and share code on GitHub.

# (Optional) Install Docker

Docker makes your software stack portable and reproducible.

sudo apt install docker.io -y
sudo usermod -aG docker $USER

docker.io → container runtime.
usermod -aG docker $USER → allows you to run Docker without sudo.

👉 Why now?
This is advanced but useful for sharing your pipelines with collaborators or running tools not available on Conda.

# ✅ Summary

After completing this setup, you’ll be able to:

Download sequencing data from NCBI/ENA.

Run quality control and trimming.

Align reads to genomes.

Process alignment (BAM/VCF) files.

Run RNA-seq differential expression analysis.

Document and share your work reproducibly.
