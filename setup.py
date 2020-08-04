from setuptools import setup, find_namespace_packages

setup(name='m3h3',
      version='0.1',
      description='An open source framework for cardiac modelling using FEniCS.',
      # url='https://github.com/ComputationalPhysiology/m3h3',
      url = "https://github.com/MartinKHovden/m3h3",
      author='Alexandra K. Diem, Martin K. Hovden',
      author_email='alexandra@simula.no, marthvdn@online.no', 
      license='BSD3',
      packages=find_namespace_packages(),
      zip_safe=False)
