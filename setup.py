# -*- coding: utf-8 -*-
"""
Setup install script for genmechanics

"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()
    
    
from os.path import dirname, isdir, join
import re
from subprocess import CalledProcessError, check_output


tag_re = re.compile(r'\btag: %s([0-9][^,]*)\b')
version_re = re.compile('^Version: (.+)$', re.M)


def get_version():
    # Return the version if it has been injected into the file by git-archive
    version = tag_re.search('$Format:%D$')
    if version:
        return version.group(1)

    d = dirname(__file__)
    
    if isdir(join(d, '.git')):
        cmd = 'git describe --tags'
        try:
            version = check_output(cmd.split()).decode().strip()[:]
        except CalledProcessError:
            raise RuntimeError('Unable to get version number from git tags')
        if version[0]=='v':
            version = version[1:]
        # PEP 440 compatibility
        if '-' in version:
            future_version = version.split('-')[0].split('.')
            future_version[-1] = str(int(future_version[-1])+1)
            future_version = '.'.join(future_version)
            number_commits = version.split('-')[1]
            version = '{}.dev{}'.format(future_version, number_commits)
            return version

    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)
            
    # Writing to file
    with open('genmechanics/version.py', 'w') as vf:
        vf.write("# -*- coding: utf-8 -*-\nversion = '{}'".format(version))
                 
    return version



setup(name='genmechanics',
#      use_scm_version={'write_to':'genmechanics/version.py'},
#      setup_requires=['setuptools_scm'],
      version = get_version(),
      description='General mechanics solver for python',
      long_description=readme(),
      keywords='General mechanics',
      url='https://github.com/Dessia-tech/genmechanics',
      author='Steven Masfaraud',
      author_email='root@dessia.tech',
      license='Creative Commons Attribution-Share Alike license',
      packages=['genmechanics'],
      package_dir={'genmechanics': 'genmechanics'},
      install_requires=['numpy', 'matplotlib', 'networkx', 'scipy', 'volmdlr'],
      classifiers=['Topic :: Scientific/Engineering','Development Status :: 3 - Alpha'])
