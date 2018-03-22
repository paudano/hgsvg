unset PYTHONPATH
virtualenv environments/python2.7
source ./environments/python2.7/bin/activate

pip install numpy
pip install pandas
pip install pysam
pip install networkx
pip install datetime
pip install argparse
pip install intervaltree
pip install biopython
pip install h5py
