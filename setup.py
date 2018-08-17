# -*- coding: utf-8 -*-
"""
Setup install script for genmechanics

"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='genmechanics',
      use_scm_version={'write_to':'genmechanics/version.py'},
      setup_requires=['setuptools_scm'],
      description='General mechanics solver for python',
      long_description=readme(),
      keywords='General mechanics',
      url='http://github.com/masfaraud/genmechanics',
      author='Steven Masfaraud',
      author_email='steven@masfaraud.fr',
      license='Creative Commons Attribution-Share Alike license',
      packages=['genmechanics'],
      package_dir={'genmechanics':'genmechanics'},
      install_requires=['numpy','matplotlib','networkx','scipy'],
      classifiers=['Topic :: Scientific/Engineering','Development Status :: 3 - Alpha'])
