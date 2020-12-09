#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

import netsyn

setup(
    name = 'netsyn',
    version = netsyn.__version__,
    packages = find_packages(),

    author = 'CEA/DRF/JACOB/GENOSCOPE/LABGeM',
    author_email = 'labgem@genoscope.cns.fr',

    description = 'NetSyn is a tool to detect conserved genomic contexts (i.e. synteny conservation) among a list of protein targets.',
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

    url = 'http://github.com/labgem/netsyn',

    classifiers = [
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],

    entry_points = {
        'console_scripts': [
            'netsyn=netsyn.netsyn:main',
            'netsyn_getINSDCFiles=netsyn.netsyn_getINSDCFiles:main',
            'netsyn_parseINSDCFiles_GetTaxonomy=netsyn.netsyn_parseINSDCFiles_GetTaxonomy:main',
            'netsyn_clusteringIntoFamilies=netsyn.netsyn_clusteringIntoFamilies:main',
            'netsyn_syntenyFinder=netsyn.netsyn_syntenyFinder:main',
            'netsyn_dataExport=netsyn.netsyn_dataExport:main'
        ]
    }
)