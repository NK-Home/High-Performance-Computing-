To run industry-standard bioinformatics tools (like samtools or bedtools) on Windows 
you must use WSL (Windows Subsystem for Linux). This creates a Linux environment inside your Windows computer, which is necessary because these tools are rarely built for Windows.



# Step 1: Install WSL and Mamba (The Speed Boost)

Bioinformatics tools have complex "dependencies" (requirements). 

Standard Anaconda is often too slow to solve them, so we use Mamba, which is a high-speed version of Conda.

- Open PowerShell: Click the Windows Start menu, type PowerShell, right-click it, and select Run as Administrator.

- Enable Linux: Type wsl --install and press Enter. Once finished, Restart your computer.

- Open Ubuntu: After restarting, search for Ubuntu in your Start menu and open it. Follow the prompts to create a simple username and password.

- Install Mamba inside Ubuntu: In the Ubuntu window, paste these two commands (one at a time):

      bash
      wget github.com
      bash Mambaforge-Linux-x86_64.sh
      Use code with caution.

- Follow the prompts: Press Enter to read the license, type yes to accept, and type yes when asked to "initialize."

- Close and reopen Ubuntu to finish - y

# Step 2: Configure "Bioconda"

Bioconda is the specialized "app store" for bioinformatics. You must tell Mamba to look there for tools.

- In the Ubuntu terminal, enter these four commands one by one:

      bash
      mamba config --add channels defaults
      mamba config --add channels bioconda
      mamba config --add channels conda-forge
      mamba config --set channel_priority strict


# Step 3: Create Your Bioinformatics Environment

This creates a "container" where your specific versions of Python and your bash tools (like samtools) will live safely without messing up your computer.

- In the Ubuntu terminal, enter:
  
      bash
      mamba create -n bio_work python=3.10 ipykernel samtools bedtools

- bio_work: The name of your environment.
- ipykernel: Essential for VS Code to "see" this environment.
- samtools & bedtools: Your first two bioinformatics bash tools.
- Type y when prompted to finish.
  
# Step 4: Link VS Code to your Linux Environment

This is the most important step for Windows users. You need to tell VS Code to "look inside" the Linux subsystem.

- Open VS Code on Windows.

- Click the Extensions icon on the left (looks like 4 squares) and install the WSL Extension.
- Click the Blue "><" icon in the very bottom-left corner of VS Code.
- Select Connect to WSL. VS Code will reload; you are now working "inside" Linux.
- Open your project folder via File > Open Folder.
- To set the Kernel: Open a .py or .ipynb file. Press Ctrl+Shift+P, type Python: Select Interpreter, and select 'bio_work': conda.
  
# Step 5: Using Bash Tools in VS Code

Because you are connected via WSL, you can use Linux bash tools directly inside VS Code.
The Integrated Terminal: Click Terminal > New Terminal. It will automatically be an Ubuntu terminal. 
Type samtools --version to see it work.

Inside a Notebook: To run a bash command inside a Python cell, start the line with an exclamation mark:

      python
      !samtools --version


# Summary for Windows Users


- Initial Setup	      The Ubuntu App (the black terminal window).
- Connecting VS Code	Click the Blue "><" button in the bottom-left of VS Code.
- Choosing the Kernel	The Command Palette (Ctrl+Shift+P) inside VS Code.
- Running Bash Tools	The Integrated Terminal at the bottom of VS Code.
