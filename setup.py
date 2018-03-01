# -*- coding: utf-8 -*-
"""
Setup install script for genmechanics

"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

def version_scheme(version):
    return '.'.join([str(i) for i in version.tag._key[1]])

def local_scheme(version):
    return ''

setup(name='genmechanics',
      use_scm_version={'version_scheme':version_scheme,'local_scheme':local_scheme},
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
