[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/TeamOverkill/MD_lab/HEAD) [![forthebadge](https://forthebadge.com/images/badges/made-with-c-plus-plus.svg)](https://forthebadge.com) [![forthebadge](https://forthebadge.com/images/badges/made-with-python.svg)](https://forthebadge.com) [![forthebadge](https://forthebadge.com/images/badges/you-didnt-ask-for-this.svg)](https://forthebadge.com)
# Laboratory Exercise 1: Molecular Dynamics
This repository contains two laboratory exercises in relationship to the topic of molecular dynamics as a part of the course _Statistical Thermodynamics and Molecular Simulation (KEMM38)_ 2021.

## Usage
To open the Notebooks, we consider two options. 1. The usage of Binder (recommended). 2. Running the Notebook on your local computer.
### 1. Binder
To open the Notebooks in Binder, simply click the Binder shield/URL in the top of the document or by clicking [here](https://mybinder.org/v2/gh/TeamOverkill/MD_lab/HEAD).
### 2. Local computer
To open the Notebooks, install Python3 via [Miniconda](https://conda.io/miniconda.html) or [Anaconda](https://www.anaconda.com/distribution/), and make sure all required packages are loaded using the following terminal commands
```bash
	conda env create -f environment.yml
	source activate MD_lab 
	jupyter-notebook
```
This will open the Jupyter Notebook Folder with the root being the folder you executed the previous commands. In Jupyter open a new terminal and execute the command `./postBuild` installing the molecular dynamics software and enabling Jupyter widgets.

## Layout
- `Part1_BarrierCrossing.ipynb` Jupyter Notebook to perform 1D molecular dynamics simulations using Python and generate plots.
- `Part2_Diffusion.ipynb` Jupyter Notebook to perform 3D molecular dynamics simulations using in-house software and generate plots.
- `include/` Folder containing premade data, figures and images.
- `MolecularOverkillEngine/` In-house developed molecular dynamics software.
- `figures/` Folder containing figures generated from the Jupyter Notebooks.
- `data` Folder containing start configurations, force field files and software input.

## Evaluation
Each student is expected to hand in a written report based on the two labs with the grades passed or not-passed. The written report is to be handed in no later than Fri 24/05 23:59 with the posibility of handing in a first draft no later than Mon 13/05 23:59 to recive comments no later than Fri 17/05.


If you have any questions, you can contact the lab responsibles on the following electronic addresses<br/>
stefan.hervo_hansen@teokem.lu.se<br/>
vidar.aspelin@teokem.lu.se<br/>
samuel.stenberg@teokem.lu.se<br/>

