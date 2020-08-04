#!/bin/bash

# install packages

# fenics
echo installing fenics...
conda install -c conda-forge fenics
echo

# m3h3
echo installing m3h3...
git clone https://github.com/MartinKHovden/m3h3.git
cd m3h3
python setup.py develop
cd ..
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
pip install h5py --no-binary=h5py
