#!/bin/bash

# Set proxy
export X509_USER_PROXY=$PWD/GRID_PROXY

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/usr/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/usr/etc/profile.d/conda.sh" ]; then
        . "/usr/etc/profile.d/conda.sh"
    else
        export PATH="/usr/bin:$PATH"
    fi
fi
unset __conda_setup

echo "[wrapper] hostname  = " `hostname`
echo "[wrapper] date      = " `date`
echo "[wrapper] linux timestamp = " `date +%s`

echo "[wrapper] ls-ing files (before running)"
ls -altrh

xrdcp XRD_CONDA_TARFILE .
xrdcp XRD_ANALYSIS_TARFILE .

mkdir higgs-dna
cd higgs-dna
mv ../higgs-dna.tar.gz .
tar -xzf higgs-dna.tar.gz
cd ..

echo "[wrapper] ls-ing files (after untarring)"
ls -altrh
ls -althr higgs-dna/

# activate env
source higgs-dna/bin/activate
export PYTHONPATH=`pwd`:$PYTHONPATH
export PATH=`pwd`/higgs-dna/bin:$PATH

echo $PATH
echo $PYTHONPATH

# Untar analysis environment
tar xf higgs_dna.tar.gz
pip install -e .

python PYTHON_FILE
