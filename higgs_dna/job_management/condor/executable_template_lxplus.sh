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

# Untar
mkdir higgs-dna
cd higgs-dna
mv ../higgs-dna.tar.gz .
tar xf higgs-dna.tar.gz
cd ..

echo "[wrapper] ls-ing files (after untarring)"
ls -altrh
ls -althr higgs-dna/

# activate env
echo "[wrapper] Activating higgs-dna env"
source higgs-dna/bin/activate
export PYTHONPATH=`pwd`:$PYTHONPATH
export PATH=`pwd`/higgs-dna/bin:$PATH

export PYTHONPATH=`pwd`:$PYTHONPATH
export PATH=`pwd`/higgs-dna/bin:$PATH

echo $PATH
echo $PYTHONPATH

which python
which python3
which pip
which pip3
python -V
python3 -V

pushd HIGGS_DNA_BASE
#pip install --upgrade pip
#pip install --upgrade setuptools
pip install -e .
popd

cd OUTPUT_DIR

python PYTHON_FILE
