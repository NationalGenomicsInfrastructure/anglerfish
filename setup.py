#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

setup(
    name='anglerfish',
    version='0.3.0',
    description='Anglerfish, a tool to demultiplex Illumina libraries from ONT data',
    author='Remi-Andre Olsen',
    author_email='remi-andre.olsen@scilifelab.se',
    url='https://github.com/remiolsen/anglerfish',
    license='MIT',
    packages = find_packages(),
    install_requires=[
        'python-levenshtein',
        'biopython'
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
