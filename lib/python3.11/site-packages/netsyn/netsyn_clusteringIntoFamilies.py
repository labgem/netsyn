#!/usr/bin/env python3

##########
# Import #
##########
from netsyn import common

import re
import sys
import os
import shutil
import subprocess
import logging
import argparse
#############
# Functions #
#############


def removeUselessFiles(dataDirectoryProcess):
    for fileName in os.listdir(dataDirectoryProcess):
        if not (fileName.endswith('.tsv') or fileName.endswith('.log')):
            os.remove(os.path.join(dataDirectoryProcess, fileName))


def mmseqs_createdb(dataDirectoryProcess, multiFasta):
    ''' create database using the mmseqs software
    '''
    logger = logging.getLogger('{}.{}'.format(
        mmseqs_createdb.__module__, mmseqs_createdb.__name__))
    dataBase_path = os.path.join(dataDirectoryProcess, 'dataBase.DB')
    logFile = os.path.join(dataDirectoryProcess, 'mmseqs_createdb.log')
    with open(logFile, 'w') as file:
        try:
            subprocess.run(['mmseqs', 'createdb', multiFasta,
                            dataBase_path], stdout=file, stderr=file, check=True)
            logger.info('Createdb completed')
        except:
            logger.error(
                'mmseqs createdb failed. Please, check {}'.format(logFile))
            exit(1)


def mmseqs_clustering(dataDirectoryProcess, params):
    ''' cluster sequences using the mmseqs software
    '''
    logger = logging.getLogger('{}.{}'.format(
        mmseqs_clustering.__module__, mmseqs_clustering.__name__))
    mmseqsTMPdirectory = os.path.join(dataDirectoryProcess, 'MMseqsTMP')
    os.mkdir(mmseqsTMPdirectory)
    dataBase_path = os.path.join(dataDirectoryProcess, 'dataBase.DB')
    cluster_path = os.path.join(dataDirectoryProcess, 'cluster.cluster')
    clustering_settings = ['mmseqs', 'cluster', dataBase_path,
                           cluster_path, mmseqsTMPdirectory,
                           '--min-seq-id', params['min_id'],
                           '-c', params['min_coverage'],
                           ]
    advanced_settings = common.readYamlAdvancedSettingsFile(
        params['advancedSettings'], common.getMMseqsDefaultSettings())
    settings_separator = '_'
    for settings in advanced_settings.values():
        for name, value in settings.items():
            clustering_settings.append(
                '--{}'.format(settings_separator.join(name.split(settings_separator)[1:])))
            clustering_settings.append(str(value))
    logFile = os.path.join(dataDirectoryProcess, 'mmseqs_clustering.log')
    with open(logFile, 'w') as file:
        try:
            subprocess.run(clustering_settings, stdout=file,
                           stderr=file, check=True)
            logger.info('Clustering completed')
        except:
            logger.error(
                'mmseqs cluster failed. Please, check {}'.format(logFile))
            exit(1)


def mmseqs_createTSV(dataDirectoryProcess, outputTSV_path):
    ''' execute the mmseqs command line 'mmseqs createtsv'
    '''
    logger = logging.getLogger('{}.{}'.format(
        mmseqs_createTSV.__module__, mmseqs_createTSV.__name__))
    dataBase_path = os.path.join(dataDirectoryProcess, 'dataBase.DB')
    cluster_path = os.path.join(dataDirectoryProcess, 'cluster.cluster')
    logFile = os.path.join(dataDirectoryProcess, 'mmseqs_createtsv.log')
    with open(logFile, 'w') as file:
        try:
            subprocess.run(['mmseqs', 'createtsv', dataBase_path,
                            dataBase_path, cluster_path, outputTSV_path
                            ], stdout=file, stderr=file, check=True)
            logger.info('CreateTSV completed')
        except:
            logger.error(
                'mmseqs createtsv failed. Please, check {}'.format(logFile))
            exit(1)


def regroup_families(tsv_file, prots_info):
    ''' create a dictionary to store families obtained by MMseqs2
    '''
    INC_FAMILY = 0
    centroid = None
    tmp_dict = {}
    lines = common.read_file(tsv_file)
    for aline in lines:
        aline = aline.strip().split('\t')
        if centroid != aline[0]:
            centroid = aline[0]
            INC_FAMILY += 1
        cds = aline[1]
        tmp_dict[cds] = INC_FAMILY

    for idx, _ in enumerate(prots_info):
        prots_info[idx]['family'] = tmp_dict[prots_info[idx]['id']
                                             ] if prots_info[idx]['id'] in tmp_dict else common.global_dict['defaultValue']
    return prots_info


def run(FASTA_FILE, PROTEINS, IDENTITY, COVERAGE, ADVANCEDSETTINGSFILENAME):
    ''' main script to run the second box of NetSyn2
    '''
    # Constants
    # common.constantsInitialiszation(args.ProjectName, args.InputFile) # depends on how the function is launched (by hand or via netsyn)
    boxName = common.global_dict['boxName']['ClusteringIntoFamilies']
    dataDirectoryProcess = os.path.join(
        common.global_dict['dataDirectory'], boxName)
    # Outputs
    families = common.global_dict['files'][boxName]['families']
    proteins_2 = common.global_dict['files'][boxName]['proteins_2']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    reportingMessages = []
    print('')
    logger.info('{} running...'.format(boxName))
    # Process
    params = {
        'min_id': str(IDENTITY),
        'min_coverage': str(COVERAGE),
        'advancedSettings': ADVANCEDSETTINGSFILENAME
    }

    if not os.path.isdir(dataDirectoryProcess):
        os.mkdir(dataDirectoryProcess)
    else:
        if 'MMseqsTMP' in os.listdir(dataDirectoryProcess):
            shutil.rmtree(os.path.join(dataDirectoryProcess, 'MMseqsTMP/'))
        removeUselessFiles(dataDirectoryProcess)

    mmseqs_createdb(dataDirectoryProcess, FASTA_FILE)
    mmseqs_clustering(dataDirectoryProcess, params)
    mmseqs_createTSV(dataDirectoryProcess, families)
    shutil.rmtree(os.path.join(dataDirectoryProcess, 'MMseqsTMP/'))
    removeUselessFiles(dataDirectoryProcess)

    prots_info = common.readJSON(
        PROTEINS, common.getProteinsParsingStepSchema())
    prots_info = regroup_families(families, prots_info)
    common.write_json(prots_info, proteins_2)

    logger.info('{} completed!'.format(boxName))

    countEachFamily = {}
    for protein in prots_info:
        if protein['family'] not in countEachFamily:
            countEachFamily[protein['family']] = 1
        else:
            countEachFamily[protein['family']] += 1

    countSingeton = len(
        [1 for _, count in countEachFamily.items() if count == 1])
    reportingMessages.append('Genomic context size processed: {}'.format(
        common.global_dict['maxGCSize']))
    reportingMessages.append(
        'Proteins number processed: {}'.format(len(prots_info)))
    reportingMessages.append('Proteins families number: {}'.format(
        len(countEachFamily)-countSingeton
    ))
    reportingMessages.append(
        'Proteins singleton number: {}'.format(countSingeton))
    reportingMessages.append(
        'Proteins cluster number: {}'.format(len(countEachFamily)))
    common.reportingFormat(logger, boxName, reportingMessages)


def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                     usage='''ClusteringIntoFamilies.py -f <fastaFileName> -ip <proteinsFile> -o <outputName> -id <ident> -mc <MinimalCoverage>''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-f', '--FastaFile', type=str,
                        required=True, help='Fasta file obtained during the previous process ParseINSDCFiles_GetTaxonomy with all proteins sequences')
    group1.add_argument('-ip', '--inputProteins', type=str,
                        required=True, help='Json file obtained during the previous process ParseINSDCFiles_GetTaxonomy containing proteins information')
    group1.add_argument('-o', '--outputName', type=str,
                        required=True, help='Output name files')
    group1.add_argument('-id', '--Identity', type=float,
                        default=0.3, help='Sequence identity.\nDefault value: 0.3')
    group1.add_argument('-cov', '--Coverage', type=float,
                        default=0.8, help='Minimal coverage allowed.\nDefault value: 0.8')

    group3 = parser.add_argument_group('Advanced settings')
    group3.add_argument('-asm', '--MMseqsAdvancedSettings', type=str,
                        help='YAML file with the advanced clustering settings to determine protein families. Settings of MMseqs2 software')

    group2 = parser.add_argument_group('logger')
    group2.add_argument('--logLevel',
                        type=str,
                        nargs='?',
                        default='INFO',
                        help='log level',
                        choices=['ERROR', 'error', 'WARNING',
                                 'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                        required=False)
    group2.add_argument('--logFile',
                        type=str,
                        nargs='?',
                        help='log file (use the stderr by default)',
                        required=False)
    return parser.parse_args(), parser


def main():
    ######################
    # Parse command line #
    ######################
    args, parser = argumentsParser()
    for key, value in {'Identity': args.Identity, 'Coverage': args.Coverage}.items():
        if not (0 < value <= 1):
            parser.error(
                'ValueError: value of --{} option must be a float number according to the condition: 0 < value <= 1'.format(key))
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
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault(
        'proteins_2', '{}_proteins_familiesStep.json'.format(args.outputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault(
        'families', '{}_families.tsv'.format(os.path.join(common.global_dict['dataDirectory'], boxName, args.outputName)))
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault(
        'report', '{}_{}_report.txt'.format(args.outputName, boxName))
    #######
    # Run #
    #######
    run(args.FastaFile, args.inputProteins, args.Identity,
        args.Coverage, args.MMseqsAdvancedSettings)


if __name__ == '__main__':
    main()
