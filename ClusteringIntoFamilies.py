#!/usr/bin/env python3

##########
# Import #
##########
import common
import re
import sys
import os
import shutil
import subprocess
import logging
#############
# Functions #
#############

def mmseqs_createdb(dataDirectoryProcess, multiFasta, prefix):
    ''' create database using the mmseqs software
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_createdb.__module__, mmseqs_createdb.__name__))
    with open('{}/{}'.format(dataDirectoryProcess, 'mmseqs_createdb.log'), 'w') as file:
        db_creation = subprocess.run(['mmseqs', 'createdb', multiFasta, '{}/{}.{}'.format(dataDirectoryProcess, prefix, 'DB')], stdout=file, stderr=file, check=True)
        logger.info('createdb - exit code: {}'.format(db_creation.returncode))
    return 0

def mmseqs_clustering(dataDirectoryProcess, params):
    ''' cluster sequences using the mmseqs software
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_clustering.__module__, mmseqs_clustering.__name__))
    mmseqsTMPdirectory = '{}/{}'.format(dataDirectoryProcess, 'MMseqsTMP')
    if os.path.isdir(mmseqsTMPdirectory):
        shutil.rmtree(mmseqsTMPdirectory)
    os.mkdir(mmseqsTMPdirectory)
    prefix = params['prefix']
    dataBase_name = '.'.join([prefix, 'DB'])
    dataBase_path = os.path.join(dataDirectoryProcess, dataBase_name)
    suffixCluster = 'cluster'
    outputCluster_name = '.'.join([prefix, 'cluster'])
    outputCluster_path = os.path.join(dataDirectoryProcess, outputCluster_name)
    with open('{}/{}'.format(dataDirectoryProcess, 'mmseqs_clustering.log'), 'w') as file:
        clust_creation = subprocess.run(['mmseqs', 'cluster', dataBase_path,
                                         outputCluster_path, mmseqsTMPdirectory,
                                         '--min-seq-id', params['min_id'],
                                         '--cov-mode', params['cov_mode'],
                                         '-c', params['min_coverage'],
                                         '--cluster-mode', str(2),
                                         '--kmer-per-seq', str(80),
                                         '--max-seqs', str(300)
                                         ], stdout=file, stderr=file, check=True)
        logger.info('clustering - exit code: {}'.format(clust_creation.returncode))
    return 0

def mmseqs_createTSV(dataDirectoryProcess, prefix):
    ''' execute the mmseqs command line 'mmseqs createtsv'
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_createTSV.__module__, mmseqs_createTSV.__name__))
    inputDB_name = '.'.join([prefix, 'DB'])
    inputCluster_name = '.'.join([prefix, 'cluster'])
    outputTSV_name = '.'.join([prefix, 'tsv'])
    inputDB_path = os.path.join(dataDirectoryProcess, inputDB_name)
    inputCluster_path = os.path.join(dataDirectoryProcess, inputCluster_name)
    outputTSV_path = os.path.join(dataDirectoryProcess, outputTSV_name)
    with open('{}/{}'.format(dataDirectoryProcess, 'mmseqs_createtsv.log'), 'w') as file:
        tsv_creation = subprocess.run(['mmseqs', 'createtsv', inputDB_path,
                                       inputDB_path, inputCluster_path, outputTSV_path
                                       ], stdout=file, stderr=file, check=True)
        logger.info('createTSV - exit code: {}'.format(tsv_creation.returncode))
    return 0

def mmseqs_runner(params, dataDirectoryProcess, multiFasta):
    ''' runs the mmseqs2 software on the multiFasta file 'ALL.faa'
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_runner.__module__, mmseqs_runner.__name__))
    logger.info('MMseqs2 running ...')
    
    logger.info('End of MMseqs2 running !')
    return 0

def regroup_families(tsv_file, prots_info):
    ''' create a dictionary to store families obtained by MMseqs2
    '''
    INC_FAMILY = 1
    centroid = None
    tmp_dict = {}
    lines = common.read_file(tsv_file)
    for aline in lines:
        aline = aline.strip().split('\t')
        if not centroid:
            centroid = aline[0]
        elif centroid != aline[0]:
            centroid = aline[0]
            INC_FAMILY += 1
        cds = int(aline[1])
        tmp_dict[cds] = INC_FAMILY

    for idx, _ in enumerate(prots_info):
        prots_info[idx]['family'] = tmp_dict[prots_info[idx]['id']] if tmp_dict[prots_info[idx]['id']] else None
    return prots_info

def get_organisms_idx(targets_info, orgs_2_Out, targets_2_Out):
    ''' get the organism index in the organisms list
    input: edited file containing all organisms information
    output: targets_info dictionary updated with a new field 'organism_idx'
    '''
    orgs = common.read_pickle(orgs_2_Out)
    for idx, organism in enumerate(orgs):
        for target_idx in organism['targets_idx']:
            targets_info[target_idx].update({
                    'organism_idx': idx
                    })
    common.write_pickle(targets_info, targets_2_Out)
    return targets_info

def run(FASTA_FILE, PROTEINS, IDENT, COVERAGE):
    ''' main script to run the second box of NetSyn2
    '''
    # Constants
    # common.constantsInitialiszation(args.ProjectName, args.InputFile) # depends on how the function is launched (by hand or via netsyn)
    boxName = common.global_dict['boxName']['ClusteringIntoFamilies']
    dataDirectoryProcess = '{}/{}'.format(common.global_dict['dataDirectory'], boxName)
    # Outputs
    proteins_2_Out = common.global_dict['files'][boxName]['proteins_2']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    print('')
    logger.info('{} running...'.format(boxName))
    # Process
    if not os.path.isdir(dataDirectoryProcess):
        os.mkdir(dataDirectoryProcess)

    params = {
        'min_id': str(IDENT),
        'cov_mode': str(1),
        'min_coverage': str(COVERAGE)
        }
    params['prefix'] = '.'.join(os.path.basename(FASTA_FILE).split('.')[:-1])

    # MMseq2 files removing (inclure dans netsyn lors nouvelle analyse ?)
    try:
        os.remove('{}/{}.{}'.format(dataDirectoryProcess, params['prefix'], 'cluster'))
    except:
        pass
    try:
        os.remove('{}/{}.{}'.format(dataDirectoryProcess, params['prefix'], 'cluster.index'))
    except:
        pass
    try:
        os.remove('{}/{}.{}'.format(dataDirectoryProcess, params['prefix'], 'tsv'))
    except:
        pass
    try:
        os.remove('{}/{}'.format(dataDirectoryProcess, 'mmseqs_createtsv.log'))
    except:
        pass

    mmseqs_createdb(dataDirectoryProcess, FASTA_FILE, params['prefix'])
    mmseqs_clustering(dataDirectoryProcess, params)
    mmseqs_createTSV(dataDirectoryProcess, params['prefix'])
    shutil.rmtree('{}/{}'.format(dataDirectoryProcess, 'MMseqsTMP/'))

    prots_info = common.read_pickle(PROTEINS)
    prots_info = regroup_families(os.path.join(dataDirectoryProcess, '.'.join([params["prefix"], 'tsv'])), prots_info)
    common.write_pickle(prots_info, proteins_2_Out)
    common.write_json(prots_info, '{}/{}'.format(dataDirectoryProcess, 'proteins_2.json'))
    logger.info('End of ClusteringIntoFamilies')

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                     usage='''ClusteringIntoFamilies.py -f <fastaFileName> -p <proteinsFile> -o <OutputName> -id <ident> -mc <MinimalCoverage>''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-f', '--FastaFile', type=str,
                        required=True, help='Fasta file obtained during the previous process ParseINSDCFiles_GetTaxonomy.')
    group1.add_argument('-p', '--Proteins', type=str,
                        required=True, help='Json file obtained during the previous process ParseINSDCFiles_GetTaxonomy containing proteins information.')
    group1.add_argument('-o', '--OutputName', type=str,
                        required=True, help='Output name files.')
    group1.add_argument('-id', '--Ident', type=float,
                        default=0.3, help='Sequence identity.\nDefault value: 0.3.')
    group1.add_argument('-mc', '--MinCoverage', type=float,
                        default=0.8, help='Minimal coverage allowed.\nDefault value: 0.8.')

    group2 = parser.add_argument_group('logger')
    group2.add_argument('--log_level',
                         type=str,
                         nargs='?',
                         default='INFO',
                         help='log level',
                         choices=['ERROR', 'error', 'WARNING', 'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                         required=False)
    group2.add_argument('--log_file',
                         type=str,
                         nargs='?',
                         help='log file (use the stderr by default)',
                         required=False)
    return parser.parse_args()

if __name__ == '__main__':
    import argparse
    ######################
    # Parse command line #
    ######################
    args = argumentsParser()
    ##########
    # Logger #
    ##########
    common.parametersLogger(args)
    #########################
    # Dependancies checking #
    #########################
    common.dependanciesChecking()
    #############
    # Constants #
    #############
    common.global_dict['dataDirectory'] = '.'
    boxName = common.global_dict['boxName']['ClusteringIntoFamilies']
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault('proteins_2', '{}/{}_proteins_2.pickle'.format(boxName, args.OutputName))
    #######
    # Run #
    #######
    run(args.FastaFile, args.Proteins, args.Ident, args.MinCoverage)
