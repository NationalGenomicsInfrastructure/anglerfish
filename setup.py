#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

setup(
    name='bio-anglerfish',
    version='0.4.2',
    description='Anglerfish, a tool to demultiplex Illumina libraries from ONT data',
    author='Remi-Andre Olsen',
    author_email='remi-andre.olsen@scilifelab.se',
    url='https://github.com/remiolsen/anglerfish',
    license='MIT',
    packages = find_packages(),
    install_requires=[
        'python-levenshtein==0.12.0',
        'biopython==1.70',
        'numpy==1.19.2'
    ],
    scripts=['./anglerfish.py'],
    zip_safe=False,
    classifiers=[
    	"Development Status :: 3 - Alpha",
    	"Environment :: Console",
    	"Intended Audience :: Developers",
    	"Intended Audience :: Healthcare Industry",
    	"Intended Audience :: Science/Research",
    	"License :: OSI Approved :: MIT License",
    	"Operating System :: POSIX :: Linux",
    	"Programming Language :: Python",
    	"Topic :: Scientific/Engineering :: Medical Science Apps."
	]
)
