#!/bin/bash

# Set proxy
export X509_USER_PROXY=$PWD/x509up_u134033

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

/eos/user/z/zhenxuan/higgs-dna/bin/python /eos/user/z/zhenxuan/hhwwgg/GluGluToRadionToHHTo2G4Q_M1000_2017/job_9/GluGluToRadionToHHTo2G4Q_M1000_2017_executable_job9.py
