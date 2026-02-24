# Python Basics Trial - Notes

## Set up 
- Download Python (most recent version) [link] (https://www.python.org/downloads/)
- Download VS code -> add python and jyoter notebook extension 

## Set up a file and an vitural enviroment
- set up a jyputer notebook file by =.ipynb file by using Cntrl+Shift+P o open command pallete and used prompt - create new jyputer notebook.
- Use cntrl+S to save it in file and give it a title. the file ends in .json file 

Connect to  Kernel 
- engine that excutes the code - top right - select kernel = allows you to connect to a python enviroment or jypter server
- choose a python enviroment 
-> use virtual enviroment rather than global enviroment = allows you to silo different projects.
  
 venv = used for general projects
 Windows - in the terminal 
-  You can also use `py -3 -m venv .venv`
python -m venv (title).venv
you create a new virtual environment, a prompt will be displayed in VS Code to allow you to select it for the workspace

- or though the command pallete
  1 - create enviroment search
  2 - select the venv and - choose an interpertur = python version
  
  
 conda = used for data science projects 
  - create this conda enviroment in the terminal => bash code in command prompt
    conda create -n env-01 python=3.9 scipy=0.15.0 numpy

 
