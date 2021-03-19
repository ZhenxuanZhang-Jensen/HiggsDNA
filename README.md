# HiggsDNA
Higgs to diphoton nanoAOD framework

## Installation

The installation procedure consists in the following steps:

 **1. Clone this repository**
```
git clone https://gitlab.cern.ch/HiggsDNA-project/HiggsDNA
cd HiggsDNA
```
**2. Install dependencies**

The necessary dependencies (listed in ```environment.yml```) can be installed manually, but the suggested way is to create a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) by running: 
```
conda env create -f environment.yml
```
Please note that the field ```python>=3.6``` will create an environment with the most recent stable version of Python. Change it to suite your needs (but still matching the requirement of Python>=3.6). 

 **3. Install ```higgs_dna```**
 
**Users** can install the package by simply running:
```
python setup.py install
```
(when a stable version will be available, it will be uploaded to the PyPI and it will be possible to install with just ```pip install higgs_dna``` without the need to clone the repository). 


For **developers**, the suggested way to install is:
```
pip install -e .
```
this prevents the need to run the installation step every time a change is performed.
