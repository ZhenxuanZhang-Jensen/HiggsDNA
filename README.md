# HiggsDNA  
Higgs to diphoton nanoAOD framework  
  
## Installation  
  
The installation procedure consists in the following steps:  
  
**1. Clone this repository**  
```  
git clone --recursive https://gitlab.cern.ch/HiggsDNA-project/HiggsDNA  
cd HiggsDNA  
```  
**2. Install dependencies**  
  
The necessary dependencies (listed in ```environment.yml```) can be installed manually, but the suggested way is to create a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/mana  
ge-environments.html) by running:  
```  
conda env create -f environment.yml  
```  
the conda env can become pretty large (multiple GB), so you may want to specify an installation location in the above step with
```
conda env create -f environment.yml -p <path to install conda env>
```
then activate the environment with
```
conda activate higgs-dna
```
If you specified an installation path, your activation command will instead be:
```
conda activate <path to install conda env>/higgs-dna
```

** lxplus-specific notes **
If you are running on `lxplus` you may run into permission errors, which can be fixed with:
```
chmod 755 -R /afs/cern.ch/user/<your_username_first_initial>/<your_username>/.conda
```
You may also want to increase your disk quota at [this link](https://resources.web.cern.ch/resources/Manage/EOS/Default.aspx), otherwise you may run out of space while installing your `conda` environment.

Please note that the field ```python>=3.6``` will create an environment with the most recent stable version of Python. Change it to suite your needs (but still matching the requirement of Python>=3.6).  

One additional package, `correctionlib`, must be installed via `pip`, rather than `conda`. Run
```
setup.sh
```
to install this script.
  
**3. Install ```higgs_dna```**  
  
Install with:
```  
pip install -e .  
```  

If you notice issues with the ```conda pack``` command for creating the tarball, try updating and cleaning your environment with (after running ```conda activate higgs-dna```):
```
conda env update --file environment.yml --prune
```

## Development Guidelines and Contribution
Here we summarize some general guidelines for developers and contributors.

### Printing Useful Analysis Info
We use the Python [logging facility](https://docs.python.org/3/library/logging.html) together with [rich](https://github.com/willmcgugan/rich) (for pretty printing) to provide useful analysis information. At the time of writing, two levels of information are supported: **INFO** and **DEBUG**.

All we have to do to print information from the code we are developing consists in:

- if not already present in the submodule, add the lines 
```
import logging
logger = logging.getLogger(__name__)
```
- type ```logger.info(message)``` or ```logger.debug(message)``` depending on which level we want our ```message``` to be displayed

Note: due to the way the logging package works, if the level is set to INFO only INFO messages are printed; if the level is set to DEBUG, both DEBUG and INFO messages are printed.

### Raising Detailed Exceptions

Raising detailed exceptions is important. Whenever it is possible, let's try to raise exception with detailed messages and possible solutions for the user to fix their code.

### Tests

A battery of tests based on [unittests](https://docs.python.org/3/library/unittest.html) is available. Before sending a PR, it is suggested to check that the new code didn't break anything by running:

```
python -m unittest -v
``` 

Please note that this is good practice even if CI is available. It is indeed a waste of time and resources to trigger a build if there is something clearly wrong that can be spotted by simply running the above mentioned command.
