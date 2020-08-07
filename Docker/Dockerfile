FROM quay.io/fenicsproject/stable

USER root

RUN apt-get -qq update && \
    apt-get -y upgrade && \
    apt-get -y install python3-scipy mercurial && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Update pip
RUN pip install --upgrade pip

# Install m3h3

RUN git clone https://github.com/MartinKHovden/m3h3.git
RUN cd m3h3; python setup.py install; cd ..

# Install cbcbeat
RUN hg clone https://bitbucket.org/meg/cbcbeat
RUN cd cbcbeat; python setup.py install; cd ..

# Install pulse
RUN pip install git+https://github.com/finsberg/pulse.git

# Install fenics-geometry
RUN pip install git+https://github.com/ComputationalPhysiology/fenics-geometry.git
