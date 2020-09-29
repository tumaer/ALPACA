#!/bin/sh
python3.6 -m venv venv
source venv/bin/activate
pip3 install -r ../python

source ./initialize_cluster.sh
python3 ../python/scripts/run_testsuite.py ./Configuration/AerConfig_short.xml ../ --testsuite-name Testsuite_Merge_Compile

deactivate
rm -rf venv