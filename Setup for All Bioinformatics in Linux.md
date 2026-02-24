# Step 1: Install "Mamba" (The Speed Boost)

Bioinformatics tools often have complex requirements that make standard Anaconda very slow. The industry standard is to use Mamba, which is a faster version of Conda.

Open your Linux Terminal (Shortcut: Ctrl+Alt+T).

Type this code and press Enter:

bash

    conda install -n base -c conda-forge mamba

    Use code with caution.

When asked Proceed ([y]/n)?, type y and press Enter.

# Step 2: Configure "Bioconda"

Bioconda is the specialized "app store" for bioinformatics tools. You must tell Mamba where to find it.

In the same terminal, enter these four commands one by one:

bash

    mamba config --add channels defaults

    mamba config --add channels bioconda

    mamba config --add channels conda-forge

    mamba config --set channel_priority strict



# Step 3: Create Your Bioinformatics Environment

Instead of installing tools one by one, you will create a "container" (environment) that holds everything you need.

In the terminal, enter:

bash

    mamba create -n bio_work python=3.10 ipykernel samtools bedtools

bio_work: The name of your environment.

ipykernel: Essential for VS Code to "see" this environment.

samtools & bedtools: Common bioinformatics bash tools.

Type y when prompted to finish the installation.

# Step 4: Link to VS Code

Now, tell VS Code to use the environment you just built.

Open VS Code.

Open your project folder via File > Open Folder.

To run Python scripts:

Press Ctrl+Shift+P to open the Command Palette.

Type Python: Select Interpreter and click it.

Select the one labeled 'bio_work': conda.

To run Jupyter Notebooks (.ipynb):

Open your notebook file.

Click Select Kernel in the top-right corner.

Select Python Environments... then choose bio_work.

# Step 5: Using Bash Tools in VS Code

Since you included samtools and bedtools, you can use them inside VS Code in two ways:

The Terminal: Click Terminal > New Terminal in the top menu. If you see (bio_work) at the start of the line, you can type samtools --version immediately.

Inside a Notebook:
To run a bash command inside a Python notebook cell, start the line with an exclamation mark:

python

    !samtools --version

Summary of Where to Enter Code
Task	Where to Enter
Setup & Installs	Your main Linux Terminal window.
Switching Tools	The Command Palette (Ctrl+Shift+P) in VS Code.
Running Scripts	The Integrated Terminal at the bottom of VS Code.
