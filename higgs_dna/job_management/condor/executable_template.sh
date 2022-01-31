function stageout {
    COPY_SRC=$1
    COPY_DEST=$2
    retries=0
    COPY_STATUS=1
    until [ $retries -ge 2 ]
    do
        echo "Stageout attempt $((retries+1)): env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 7200 --verbose --checksum ADLER32 ${COPY_SRC} ${COPY_DEST}"
        env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 7200 --verbose --checksum ADLER32 ${COPY_SRC} ${COPY_DEST}
        COPY_STATUS=$?
        if [ $COPY_STATUS -ne 0 ]; then
            echo "Failed stageout attempt $((retries+1))"
        else
            echo "Successful stageout with $retries retries"
            break
        fi
        retries=$[$retries+1]
        echo "Sleeping for 30s"
        sleep 30s
    done
    if [ $COPY_STATUS -ne 0 ]; then
        echo "Removing output file because gfal-copy crashed with code $COPY_STATUS"
        env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-rm --verbose ${COPY_DEST}
        REMOVE_STATUS=$?
        if [ $REMOVE_STATUS -ne 0 ]; then
            echo "Uhh, gfal-copy crashed and then the gfal-rm also crashed with code $REMOVE_STATUS"
            echo "You probably have a corrupt file sitting on hadoop now."
            exit 1
        fi
    fi
}

echo "[wrapper] hostname  = " `hostname`
echo "[wrapper] date      = " `date`
echo "[wrapper] linux timestamp = " `date +%s`

echo "[wrapper] ls-ing files (before running)"
ls -altrh

# Untar
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

which python
which python3
which pip
which pip3
python -V
python3 -V

# Untar analysis environment
tar xf higgs_dna.tar.gz
pip install -e .

# Run analysis script
python PYTHON_FILE

echo "[wrapper] ls-ing files (after running)"
ls -altrh

for x in `ls *.parquet`
do
    stageout "file://`pwd`/$x" BATCH_OUTPUT_DIR/"$x"
done

stageout "file://`pwd`/SUMMARY_FILE" BATCH_OUTPUT_DIR/SUMMARY_FILE
