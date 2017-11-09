from distutils.core import setup

from setuptools import find_packages

from mviri.mviri_tools import tools_plotting
from mviri.mviri_l10 import read_imag2tg

setup(name='mviri_lib',
      version='1.0.0',
      description='EUMETSAT python tools for reading level 1.0 and 1.5 data',
      author='Frank Ruethrich',
      author_email='frank.ruethrich@eumetsat.int',
      url='http://www.eumetsat.int',
      packages=find_packages(),
      install_requires=['numpy','matplotlib'])