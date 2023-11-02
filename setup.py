#!/usr/bin/env python
"""
Anglerfish is a tool designed to demultiplex Illumina libraries sequenced on Oxford Nanopore flowcells. The primary purpose for this would be to do QC, i.e. to check pool balancing, assess contamination, library insert sizes and so on.

Install with pip:

    pip install bio-anglerfish

Or from bioconda:

    conda install -c bioconda anglerfish
"""
from setuptools import setup, find_packages
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

version='0.6.0'

setup(
    name='bio-anglerfish',
    version=version,
    description='Anglerfish, a tool to demultiplex Illumina libraries from ONT data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Remi-Andre Olsen',
    author_email='remi-andre.olsen@scilifelab.se',
    url='https://github.com/remiolsen/anglerfish',
    license='MIT',
    python_requires=">=3.7",
    packages = find_packages(),
    package_data = {"":["config/adaptors.yaml"]},
    install_requires=[
        'python-levenshtein==0.23.0',
        'biopython==1.79',
        'numpy==1.22.0',
        'pyyaml==6.0'
    ],
    entry_points={
        "console_scripts": [
            "anglerfish=anglerfish.anglerfish:anglerfish",
        ],
    },
    zip_safe=False,
    classifiers=[
    	"Development Status :: 5 - Production/Stable",
    	"Environment :: Console",
    	"Intended Audience :: Developers",
    	"Intended Audience :: Healthcare Industry",
    	"Intended Audience :: Science/Research",
    	"License :: OSI Approved :: MIT License",
    	"Operating System :: POSIX :: Linux",
    	"Programming Language :: Python",
        "Topic :: Scientific/Engineering",
    	"Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
