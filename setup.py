#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

import netsyn2

setup(
    name = 'netsyn',
    version = netsyn2.__version__,
    packages = find_packages(),

    author = 'CEA/JACOB/GENOSCOPE/LABGeM',
    author_email = '???@genoscope.cns.fr',

    description = 'NetSyn2 is a tool implemented to detect conserved syntenies around listed targets.',
    long_description = open('README.md').read(),

    install_requires = [
        'pyyaml',
        'python-igraph',
        'jsonschema',
        'networkx',
        'markov_clustering',
        'matplotlib',
        'urllib3',
        'biopython',
        'scipy==1.1.0'
    ],

    include_package_data = True,

    url = 'http://github.com/labgem/???',

    classifiers = [
        'Programming Language :: Python :: 3.7',
        'License :: ???',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],

    entry_points = {
        'console_scripts': [
            'netsyn=netsyn2.netsyn:main',
            'netsyn_getINSDCFiles=netsyn2.netsyn_getINSDCFiles:main',
            'netsyn_parseINSDCFiles_GetTaxonomy=netsyn2.netsyn_parseINSDCFiles_GetTaxonomy:main',
            'netsyn_clusteringIntoFamilies=netsyn2.netsyn_clusteringIntoFamilies:main',
            'netsyn_syntenyFinder=netsyn2.netsyn_syntenyFinder:main',
            'netsyn_dataExport=netsyn2.netsyn_dataExport:main'
        ]
    }
)