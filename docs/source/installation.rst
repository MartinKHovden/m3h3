*****************
Installation
*****************
There are different ways to install m3h3 depending on if you want to use the 
framework as a user or a developer. The first part consentrates on how to install 
m3h3 as a user. The second part is about how to install as a developer. 

=====================================
Installation using Docker containers
=====================================
The M3H3 package can be installed using Docker. This will set up a container 
with dolfin, cbcbeat, fenics, pulse, geometry and m3h3. Run all commands in the 
terminal. 

#. Download docker by following the instructions on https://www.docker.com/get-started.

#. Clone the `m3h3 repository <https://github.com/MartinKHovden/m3h3>`_ to your system 
   and move into the m3h3 main folder. 

#. Build the docker image using the Dockerfile by running

    .. code-block:: python 

        docker build -t m3h3 docker/

#. Run the docker container from the m3h3 folder to share folder with the container.
    
    .. code-block:: python 

        docker run -ti -v $(pwd):/home/fenics/shared --name m3h3_container m3h3

#. Next install m3h3 in the container 

    .. code-block:: python 

        python3  setup.py install 

#. Now you have m3h3 installed as well as dolfin, fenics, cbcbeat, pulse and geometry. Enjoy!

#. Now you can exit the container using the exit command 

    .. code-block:: python 

        exit 

   This will send you back to your original system. 

#. To restart the container, use the start and exec commands

    .. code-block:: python 

        docker start m3h3 
        docker exec -ti -u fenics m3h3 /bin/bash -l
   
   This will send you back into the docker shell. 
    
    
======================================
Installation using Anaconda
======================================
The M3H3 package can also be installed using Anaconda. This requires 
a linux system. If you are working on Windows, set up `Windows Subsystem 
for Linux <https://docs.microsoft.com/en-us/windows/wsl/wsl2-index>`_, and run the same steps as below. 

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

The first steps are similar to what was done for the regular users. The main 
difference is that you install m3h3 in developer mode instead of 
doing a regular install. 

#. Download docker by following the instructions on https://www.docker.com/get-started. 

#. Clone the `m3h3 repository <https://github.com/MartinKHovden/m3h3>`_ to your system 
   and move into the m3h3 main folder. 

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
   You can then change the files in the m3h3 folder and the changes will immediately 
   take action in the terminal within the docker container. 

#. The container can be exited by using the exit command  

    .. code-block:: python 

        exit 

   This will send you back to your original system. 

#. To restart the container, use the start and exec commands 

    .. code-block:: python 

        docker start m3h3 
        docker exec -ti -u fenics m3h3 /bin/bash -l
   
   This will send you back into the docker shell and development mode. 


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