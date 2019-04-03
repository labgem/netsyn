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
    with open(os.path.join(dataDirectoryProcess, 'mmseqs_createdb.log'), 'w') as file:
        db_creation = subprocess.run(['mmseqs', 'createdb', multiFasta, '{}.{}'.format(os.path.join(dataDirectoryProcess, prefix), 'DB')], stdout=file, stderr=file, check=True)
        logger.info('createdb - exit code: {}'.format(db_creation.returncode))
    return 0

def mmseqs_clustering(dataDirectoryProcess, params):
    ''' cluster sequences using the mmseqs software
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_clustering.__module__, mmseqs_clustering.__name__))
    mmseqsTMPdirectory = os.path.join(dataDirectoryProcess, 'MMseqsTMP')
    if os.path.isdir(mmseqsTMPdirectory):
        shutil.rmtree(mmseqsTMPdirectory)
    os.mkdir(mmseqsTMPdirectory)
    prefix = params['prefix']
    dataBase_name = '.'.join([prefix, 'DB'])
    dataBase_path = os.path.join(dataDirectoryProcess, dataBase_name)
    suffixCluster = 'cluster'
    outputCluster_name = '.'.join([prefix, 'cluster'])
    outputCluster_path = os.path.join(dataDirectoryProcess, outputCluster_name)
    with open(os.path.join(dataDirectoryProcess, 'mmseqs_clustering.log'), 'w') as file:
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
    with open(os.path.join(dataDirectoryProcess, 'mmseqs_createtsv.log'), 'w') as file:
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
        tmp_dict[str(cds)] = INC_FAMILY

    for idx, _ in enumerate(prots_info):
        prots_info[idx]['family'] = tmp_dict[prots_info[idx]['id']] if tmp_dict[prots_info[idx]['id']] else None
    return prots_info

def run(FASTA_FILE, PROTEINS, IDENT, COVERAGE):
    ''' main script to run the second box of NetSyn2
    '''
    # Constants
    # common.constantsInitialiszation(args.ProjectName, args.InputFile) # depends on how the function is launched (by hand or via netsyn)
    boxName = common.global_dict['boxName']['ClusteringIntoFamilies']
    dataDirectoryProcess = os.path.join(common.global_dict['dataDirectory'], boxName)
    # Outputs
    proteins_2 = common.global_dict['files'][boxName]['proteins_2']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    reportingMessages = []
    print('')
    logger.info('{} running...'.format(boxName))
    # Process
    params = {
        'min_id': str(IDENT),
        'cov_mode': str(1),
        'min_coverage': str(COVERAGE)
        }
    params['prefix'] = '.'.join(os.path.basename(FASTA_FILE).split('.')[:-1])

    if not os.path.isdir(dataDirectoryProcess):
        os.mkdir(dataDirectoryProcess)
    elif 'MMseqsTMP/' in os.listdir(dataDirectoryProcess):
        shutil.rmtree(os.path.join(dataDirectoryProcess, 'MMseqsTMP/'))
        for fileName in os.listdir(dataDirectoryProcess):
            if (re.match(params['prefix'],fileName) and not fileName.endswith('.tsv')):
                os.remove(os.path.join(dataDirectoryProcess, fileName))

    mmseqs_createdb(dataDirectoryProcess, FASTA_FILE, params['prefix'])
    mmseqs_clustering(dataDirectoryProcess, params)
    mmseqs_createTSV(dataDirectoryProcess, params['prefix'])

    shutil.rmtree(os.path.join(dataDirectoryProcess, 'MMseqsTMP/'))

    prots_info = common.readJSON(PROTEINS)
    prots_info = regroup_families(os.path.join(dataDirectoryProcess, '.'.join([params["prefix"], 'tsv'])), prots_info)

    common.write_json(prots_info, proteins_2)
    for fileName in os.listdir(dataDirectoryProcess):
        if re.match(params['prefix'],fileName) and not fileName.endswith('.tsv'):
            os.remove(os.path.join(dataDirectoryProcess, fileName))
    logger.info('{} completed!'.format(boxName))
    countEachFamily = {}
    for protein in prots_info:
        if protein['family'] not in countEachFamily:
            countEachFamily[protein['family']] = 1
        else:
            countEachFamily[protein['family']] += 1
    countSingeton = len([1 for _, count in countEachFamily.items() if count == 1])
    reportingMessages.append('Genomic context size processed: {}'.format(common.global_dict['maxGCSize']))
    reportingMessages.append('Proteins number processed: {}'.format(len(prots_info)))
    reportingMessages.append('Proteins families number: {}'.format(
        len(countEachFamily)-countSingeton
    ))
    reportingMessages.append('Proteins singleton number: {}'.format(countSingeton))
    reportingMessages.append('Proteins cluster number: {}'.format(len(countEachFamily)))
    common.reportingFormat(logger, boxName, reportingMessages)

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                     usage='''ClusteringIntoFamilies.py -f <fastaFileName> -p <proteinsFile> -o <OutputName> -id <ident> -mc <MinimalCoverage>''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-f', '--FastaFile', type=str,
                        required=True, help='Fasta file obtained during the previous process ParseINSDCFiles_GetTaxonomy with all proteins sequences')
    group1.add_argument('-p', '--Proteins', type=str,
                        required=True, help='Json file obtained during the previous process ParseINSDCFiles_GetTaxonomy containing proteins information')
    group1.add_argument('-o', '--OutputName', type=str,
                        required=True, help='Output name files')
    group1.add_argument('-id', '--Ident', type=float,
                        default=0.3, help='Sequence identity.\nDefault value: 0.3')
    group1.add_argument('-mc', '--MinCoverage', type=float,
                        default=0.8, help='Minimal coverage allowed.\nDefault value: 0.8')

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
    return parser.parse_args(), parser

if __name__ == '__main__':
    import argparse
    ######################
    # Parse command line #
    ######################
    args, parser = argumentsParser()
    for key, value in {'Ident':args.Ident, 'MinCoverage':args.MinCoverage}.items():
        if not (0 < value <= 1):
            parser.error('ValueError: value of --{} option must be a float number according to the condition: 0 < value <= 1'.format(key))
    ##########
    # Logger #
    ##########
    common.parametersLogger(args)
    #########################
    # Dependancies checking #
    #########################
    common.dependenciesChecking()
    #############
    # Constants #
    #############
    common.global_dict['dataDirectory'] = '.'
    boxName = common.global_dict['boxName']['ClusteringIntoFamilies']
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault('proteins_2', '{}_proteins_familiesStep.json'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('report', '{}_{}_report.txt'.format(args.OutputName, boxName))
    #######
    # Run #
    #######
    run(args.FastaFile, args.Proteins, args.Ident, args.MinCoverage)
