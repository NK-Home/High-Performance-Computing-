1️⃣ Install VS Code on Windows

Go to https://code.visualstudio.com/
 and download the installer.

Install VS Code with default options (you can enable “Add to PATH” to make it easier later).

Why: VS Code is a lightweight but powerful editor that can connect to WSL (Windows Subsystem for Linux) to run Ubuntu commands seamlessly.

2️⃣ Install the WSL extension in VS Code

Open VS Code.

Go to the Extensions tab (or press Ctrl+Shift+X).

Search for “Remote - WSL”.

Click Install.

Why: This extension lets VS Code use your Ubuntu (WSL) environment as the backend for running code, terminals, and Jupyter.

3️⃣ Open VS Code connected to your WSL environment

Press Ctrl+Shift+P to open the command palette.

Type Remote-WSL: New Window and select it.

VS Code will open a new window connected to Ubuntu.

Open your project folder (e.g., ~/bioinformatics) inside this WSL session.

Tip: You can use the terminal in VS Code (Ctrl+` ) and it will run Ubuntu commands exactly as in your terminal.

4️⃣ Install Jupyter in your bioinfo environment

Make sure your bioinfo Conda environment is active:

conda activate bioinfo
mamba install jupyterlab


Why: This lets you run notebooks directly inside VS Code using your bioinfo environment with all your tools (FastQC, samtools, etc.).

5️⃣ Open a Jupyter Notebook in VS Code

In VS Code, press Ctrl+Shift+P → Python: Select Interpreter → choose your bioinfo Conda environment.

Create a new .ipynb notebook.

You can now write Python code, import packages like pandas or matplotlib, and run analysis on files stored in your WSL environment (~/bioinformatics/raw_data, etc.).

6️⃣ Access the terminal and run bioinformatics tools

In VS Code, open the integrated terminal (Ctrl+` ).

You can run all terminal commands you’ve been practicing:

fastqc raw_data/SRR2584863_1.fastq -o results/
multiqc results/
bwa mem reference/ecoli.fa raw_data/SRR2584863_1.fastq raw_data/SRR2584863_2.fastq | samtools view -bS - > results/ecoli.bam


Everything is live and interactive, plus you can combine it with Python analysis in the notebook.

7️⃣ Optional: Open MultiQC reports automatically

If you want to open HTML reports from WSL in Windows:

explorer.exe results/multiqc_report.html


Or if you’re on Linux desktop:

xdg-open results/multiqc_report.html

✅ Summary of Benefits

VS Code + WSL = one interface for everything

Run terminal commands, Python scripts, and Jupyter notebooks side by side.

Your files live in Ubuntu (~/bioinformatics) but you can view/edit them on Windows.

Full integration with Git, VS Code extensions, and debugging.
