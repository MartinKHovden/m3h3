#############
Installation
#############

In this section different methods for obtaining m3h3 is presented. 

=====================================
Installation using Docker containers
=====================================
The M3H3 package can be installed using Docker. This will set up a container 
with dolfin, cbcbeat, fenics, pulse, geometry and m3h3. Run all commands in the 
terminal. 

Should change docker file for development. 

#. Download docker from docker.com for your platform.

#. Clone

#. Build the docker image using the Dockerfile by running

    .. code-block:: python 

        docker build -t m3h3 docker/

#. The next steps are similar to before. Run the docker container from the m3h3 folder to share this with the container.
    
    .. code-block:: python 

        docker run -ti -v $(pwd):/home/fenics/shared --name m3h3_container m3h3

#. Next install m3h3 in the container 

    .. code-block:: python 

        python3  setup.py install 

#. Now you have m3h3 installed as well as dolfin, fenics, cbcbeat, pulse and geometry. Enjoy!

#. Now you can exit the container using 

    .. code-block:: python 

        exit 

   This will send you back to your original system. 

#. To restart the container, use

    .. code-block:: python 

        docker start m3h3 
        docker exec -ti -u fenics m3h3 /bin/bash -l
   
   which will send you back into the docker shell. 
    
    
======================================
Installation using Anaconda
======================================
The M3H3 package can also be installed using Anaconda. This requires 
a linux system. If you are working on Windows, set up windows subsystem 
for linux, and run the same steps as below. 

#. First, make sure that you have conda installed on your system. 

#. Next, set up a conda environment 

    .. code-block:: python 

        conda create --name m3h3 python=3.7.3

#. Activate the environment using 

    .. code-block:: python 

        conda activate m3h3 
    
#. Run the bash script to install the packages 

    .. code-block:: python 

        bash install.sh 

#. Now you have set up an environment with fenics, dolfin, cbcbeat, pulse 
   geometry, and m3h3. Enjoy! 

==================================
Installation for developers 
==================================
Both Docker and Anaconda can be used for developing.

Docker 
++++++++++

#. Download docker from docker.com for your platform. 

#. Build the docker image using the Dockerfile by running

    .. code-block:: python 

        docker build -t m3h3 docker/DockerfileDevelop

#. Next run the docker container from the m3h3 folder to share this with the container.
    
    .. code-block:: python 

        docker run -ti -v $(pwd):/home/fenics/shared --name m3h3_container m3h3

#. Next install m3h3 in the container 

    .. code-block:: python 

        python3 setup.py develop

#. Now you have m3h3 installed as well as dolfin, fenics, cbcbeat, pulse and geometry. Enjoy!

#. Now you can exit the container using 

    .. code-block:: python 

        exit 

   This will send you back to your original system. 

#. To restart the container, use

    .. code-block:: python 

        docker start m3h3 
        docker exec -ti -u fenics m3h3 /bin/bash -l
   
   which will send you back into the docker shell. 


Anaconda 
++++++++++++
#. First, set up an environment

    .. code-block:: python 

        conda create --name m3h3-develop python=3.7.3
    
#. Activate the environment using 

    .. code-block:: python 

        conda activate m3h3-develop 

#. Run the bash script to install the dependencies

    .. code-block:: python 

        bash install_dev.sh 

#. Fork the m3h3 depository to your own repo and clone it  

    .. code-block:: python 

        git clone  https://github.com/YOUR_GITHUB/m3h3.git
        cd m3h3 
        python setup.py develop 