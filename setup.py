#!/usr/bin/env python
install_requires = ['numpy','netcdf4'] 
tests_require = ['nose','coveralls']
# %%
import subprocess
from setuptools import setup,find_packages
from numpy.distutils.core import setup,Extension

ext = [Extension(name='fclimdex',
                 sources=['fclimdex.f'],
                 f2py_options=['--quiet'])]


setup(name='FClimDex',
      packages=find_packages(),
      python_requires='>=3.6',
      install_requires=install_requires,
      extras_require={'tests':tests_require,
                      'plot':['matplotlib','basemap']},
      tests_require=tests_require,
#      ext_modules=ext,
      version='0.1.0',
      author="Michael Hirsch, Ph.D.",
      url="https://github.com/scivision/geo2mag",
	  )
# %%
subprocess.check_call(['gfortran','fclimdex.f','-o','FClimdex'])
