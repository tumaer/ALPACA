#!/bin/sh
python3.6 -m venv venv
source venv/bin/activate
pip3 install --upgrade setuptools
pip3 install ../python

source ../Utilities/Clusters/LinuxCluster/modules.txt
python3 ../python/scripts/run_testsuite.py ./Configuration/.merge-compile.xml ../ --testsuite-name Testsuite_Merge_Compile

deactivate
rm -rf venv
