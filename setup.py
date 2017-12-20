#!/usr/bin/env python
install_requires = ['numpy','netcdf4','xarray','pandas']
tests_require = ['nose','coveralls']
# %%
import subprocess
from setuptools import setup,find_packages

setup(name='FClimDex',
      packages=find_packages(),
      python_requires='>=3.6',
      install_requires=install_requires,
      extras_require={'tests':tests_require,
                      'plot':['matplotlib','basemap']},
      tests_require=tests_require,
      version='0.2.0',
      author="Michael Hirsch, Ph.D.",
      url="https://github.com/scivision/fclimdex",
	  )
# %%
subprocess.check_call(['gfortran','fclimdex.f','-o','FClimdex'])
