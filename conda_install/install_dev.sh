#!/bin/bash

# install packages

# fenics
echo installing fenics...
conda install -c conda-forge fenics
echo

# pulse
echo installing dependencies
echo installing pulse...
pip install git+https://github.com/finsberg/pulse.git
echo

# fenics-geometry
echo installing fenics-geometry...
pip install git+https://github.com/ComputationalPhysiology/fenics-geometry.git
echo fenics-geometry installed successfully
echo

# cbcbeat
echo installing cbcbeat...
hg clone https://bitbucket.org/meg/cbcbeat
cd cbcbeat
python setup.py install
cd ..
echo

# h5py
echo installing h5py...
pip uninstall h5py
conda install h5py

