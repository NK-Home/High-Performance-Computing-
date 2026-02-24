## Step-by-Step Guide: Set Up Miniconda + Conda-Forge + Jupyter in VS Code

# ⭐ Step 1 — Install Miniconda

Go here:

https://docs.conda.io/projects/miniconda/en/latest/

Install the Windows version.

# ⭐ Step 2 — Open the correct terminal (important)

Use Miniconda Prompt, not:

Git Bash

PowerShell

CMD

WSL (unless doing Linux-only work)

You should see:

`(base) C:\Users\...`

# ⭐ Step 3 — Add conda-forge as your main channel

Run this:

`conda config --add channels conda-forge`

`conda config --set channel_priority strict`


This tells Conda:

👉 prefer scientific packages from conda-forge

This is EXACTLY what most bioinformaticians do.

# ⭐ Step 4 — Create a bioinformatics-friendly environment

Example:

`conda create -n bioinfo_env python=3.11`


Activate it:

`conda activate bioinfo_env`


Prompt changes to:

`(bioinfo_env) C:\Users\...>`

# ⭐ Step 5 — Install scientific Python packages + Jupyter Notebook

`conda install jupyter pandas numpy scipy seaborn matplotlib`


All from condaforge automatically.

# ⭐ Step 6 — Install the Jupyter kernel (so VS Code can see it)

`python -m ipykernel install --user --name bioinfo_env`

# ⭐ Step 7 — Open VS Code through the same terminal
`code .`


Now VS Code opens in the correct environment.

# ⭐ Step 8 — Create a Jupyter Notebook

In VS Code:

File → New File

Save as: analysis.ipynb

# ⭐ Step 9 — Select your kernel

Top right → Select Kernel → choose:

bioinfo_env

# ⭐ Step 10 — Test

`import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
print("Jupyter working in Miniconda + conda-forge!")`


If it prints → you’re fully set up.

# ⭐ Step 11 - Coming Back to the files 

when re=opening the notebook, select the correct kernel. 

Python will offload packages - create a .txt file called requirements.txt

List in the txt file, the packages you need. 

At the very beginneing of your notebook, open a cell, and run the following code to install the packages your script depeends on

`%pip install -r requirements.txt`
