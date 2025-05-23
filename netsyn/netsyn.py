#!/usr/bin/env python3

##########
# Import #
##########

from netsyn import common
from netsyn import netsyn_getINSDCFiles, netsyn_parseINSDCFiles_GetTaxonomy, netsyn_clusteringIntoFamilies, netsyn_syntenyFinder, netsyn_dataExport, EndNetSynAnalysis

import argparse
import os
import shutil
import logging
import re
import sys
import filecmp
import glob
import yaml

#############
# Functions #
#############


def checkArguments(parser):
    '''
    Arguments checking.
    '''
    args = parser.parse_args()

    if not args.UniProtACList:
        if not args.CorrespondencesFile:
            parser.error(
                'missing option: please specify the input option (-u/--UniProtACList and/or -c/--CorrespondencesFile)')
        if args.conserveDownloadedINSDC:
            parser.error(
                'ambigous option: --conserveDownloadedINSDC option requires -u/--UniProtACList option. Please specify the -u/--UniProtACList option')

    if args.GroupingOnLabel and args.GroupingOnTaxonomy:
        parser.error(
            'GroupingOnLabel and GroupingOnTaxonomy are incompatible options. Please choose one of two options')

    if (args.GroupingOnLabel or args.GroupingOnTaxonomy) and not args.ClusteringMethod:
        parser.error('The grouping clustering functionality can\'t be computed with every clustering method. Please choose a unique clustering method (-cm/--ClusteringMethod)')

    if args.ClusteringMethod and not (args.GroupingOnLabel or args.GroupingOnTaxonomy):
        parser.error('The selection of a clustering algorithm is available only when grouping clustering is enabled (see section Grouping Clustering settings: netsyn -h/--help)')

    if args.GroupingOnLabel and not args.MetaDataFile:
        parser.error('Please specify the --MetaDataFile option')

    for key, value in {'Identity': args.Identity, 'Coverage': args.Coverage}.items():
        if not (0 < value <= 1):
            parser.error(
                'ValueError: value of --{} option must be a float number according to the condition: 0 < value <= 1'.format(key))

    return args


def argumentsParser():
    '''
    Arguments parsing.
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                     usage='netsyn -u <UniProtAC.list> -o <OutputDirName> [-md <metadataFile>, --conserveDownloadedINSDC]\n\
       netsyn -u <UniProtAC.list> -c <CorrespondencesFileName> -o <OutputDirName> [-md <metadataFile>, --conserveDownloadedINSDC]\n\
       netsyn -c <CorrespondencesFileName> -o <OutputDirName> [-md <metadataFile>]',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-v', '--version', action='version',
                        version=common.global_dict['version'])

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-u', '--UniProtACList', type=str,
                        help='UniProt accession list input(cf: wiki)')
    group1.add_argument('-c', '--CorrespondencesFile', type=str,
                        help='Correspondence entry file between: protein_AC/nucleic_AC/nucleic_File_Path (cf: wiki)')
    group1.add_argument('-o', '--OutputDirName', type=str, required=True,
                        help='Output directory name, used to define the project name')
    group1.add_argument('-md', '--MetaDataFile', type=str,
                        help='File containing metadata')
    group1.add_argument('-pd', '--ProjectDescription', type=str,
                        default='No description', help='The project description')
    group1.add_argument('-np', '--newProject', action='store_true',
                        help='Force the creation of a new projet. Overwrite the project of the same name')
    group1.add_argument('--conserveDownloadedINSDC', action='store_true',
                        help='Conserve the downloaded files from GetINSDC. Require of --UniProtACList option')
    group1.add_argument('--IncludedPseudogenes', action='store_true',
                        help='CDS annotated as pseudogenes are considered as part of the genomic context')

    group4 = parser.add_argument_group(
        'Clustering settings to determine protein families')
    group4.add_argument('-id', '--Identity', type=float,
                        default=0.3, help='Sequence identity.\nDefault value: 0.3')
    group4.add_argument('-cov', '--Coverage', type=float,
                        default=0.8, help='Minimal coverage allowed.\nDefault value: 0.8')
    group4.add_argument('-mas', '--MMseqsAdvancedSettings', type=str,
                        help='YAML file with the advanced clustering settings to determine protein families. Settings of MMseqs2 software')

    group3 = parser.add_argument_group('Synteny settings')
    group3.add_argument('-ws', '--WindowSize', type=int,
                        default=common.global_dict['maxGCSize'],
                        help='Window size of genomic contexts to compare (target gene included).\nDefault value: {}.'.format(
                            common.global_dict['maxGCSize']),
                        choices=common.widowsSizePossibilities(common.global_dict['minGCSize'], common.global_dict['maxGCSize']))
    group3.add_argument('-sg', '--SyntenyGap', type=int, default=3,
                        help='Number of genes allowed between two genes in synteny.\nDefault value: 3')
    group3.add_argument('-ssc', '--SyntenyScoreCutoff', type=float,
                        default=common.global_dict['sscDefault'], help='Define the Synteny Score Cutoff to conserved.\nDefault value: >= {}'.format(common.global_dict['sscDefault']))
    group3.add_argument('-cas', '--ClusteringAdvancedSettings', type=str,
                        help='YAML file with the advanced clustering settings to determine synteny NetSyn. Settings of clusterings methods')

    group6 = parser.add_argument_group('Grouping clustering settings')
    group6.add_argument('-cm', '--ClusteringMethod', type=str,
                        choices=['MCL', 'Infomap', 'Louvain', 'WalkTrap'],
                        default=None,
                        help='Clustering method choices: MCL (small graph), Infomap (medium graph), Louvain (medium graph) or WalkTrap (big  graph).\nDefault value: MCL')
    group6.add_argument('-gl', '--GroupingOnLabel', type=str,
                        help='Label of the metadata column on which the redundancy will be computed (Incompatible with --GroupingOnTaxonomy option.)')
    group6.add_argument('-gt', '--GroupingOnTaxonomy', type=str, choices=common.global_dict['desired_ranks_lineage'].keys(),
                        help='Taxonomic rank on which the redundancy will be computed. (Incompatible with --GroupingOnLabel option)')

    group7 = parser.add_argument_group('logger')
    group7.add_argument('--logLevel',
                        type=str,
                        default='INFO',
                        help='log level',
                        choices=['ERROR', 'error', 'WARNING',
                                 'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                        required=False)
    group7.add_argument('--logFile',
                        type=str,
                        help='log file (use the stderr by default). Disable the log colors.',
                        required=False)

    args = checkArguments(parser)
    return(args)


def writeCurrentVersion(version, fileName):
    '''
    Saved the current version in a file.
    '''
    with open(fileName, 'w') as file:
        file.write(version)


def getLatestVersionUsed(fileName):
    '''
    Return the latest version used in the previous analysis.
    '''
    if os.path.isfile(fileName):
        with open(fileName, 'r') as file:
            version = file.read()
    else:
        version = common.global_dict['version']
    return version


def checkSameVersion(fileName):
    '''
    Exit if the latest version and the current version are different.
    '''
    logger = logging.getLogger(checkSameVersion.__name__)
    currentVersion = common.global_dict['version']
    latestVersion = getLatestVersionUsed(fileName)
    if not currentVersion == latestVersion:
        logger.error('The project has been launched with an anterior version (version {}). Please, create a new project or use -np option to run the project with the current version {}'.format(latestVersion, currentVersion))
        exit(1)


def orderDefinition(uniprotACList, correspondingFile):
    '''
    Setting the code execution order.
    '''
    if correspondingFile and not uniprotACList:
        order = [
            common.global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'],
            common.global_dict['boxName']['ClusteringIntoFamilies'],
            common.global_dict['boxName']['SyntenyFinder'],
            common.global_dict['boxName']['DataExport'],
            common.global_dict['boxName']['EndNetSynAnalysis']
        ]
    else:
        order = [
            common.global_dict['boxName']['GetINSDCFiles'],
            common.global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'],
            common.global_dict['boxName']['ClusteringIntoFamilies'],
            common.global_dict['boxName']['SyntenyFinder'],
            common.global_dict['boxName']['DataExport'],
            common.global_dict['boxName']['EndNetSynAnalysis']

        ]
    return order


def saveInputFiles(inputI, inputII):
    '''
    Saves the input files.
    '''
    if inputI:
        if not common.checkFilledFile(inputI):
            shutil.copyfile(inputI, common.global_dict['uniprotACListSaved'])
        else:
            exit(1)
    if inputII:
        if not common.checkFilledFile(inputII):
            shutil.copyfile(
                inputII, common.global_dict['correspondingFileSaved'])
        else:
            exit(1)


def addsAvancedSettings(yaml, step, settingsFileName, defaultValue):
    '''
    Adds the advanced settings from the file of settings or deffault settings into YAML.
    '''
    advancedSettings = common.readYamlAdvancedSettingsFile(
        settingsFileName, defaultValue)
    for settings in advancedSettings.values():
        for name, value in settings.items():
            yaml[step][name] = value


def createYamlSettingsFile(args, ORDERBOX):
    '''
    YAML settings file creation. Used for resumption of analysis.
    '''
    YAML_File = {}
    # Depending on whether execution start from GetINSDCFiles or ClusteringIntoFamilies
    gap = 0 if args.UniProtACList else 1
    with open(common.global_dict['settingsFileName'], 'w') as outfile:
        YAML_File['General informations'] = {}
        YAML_File['General informations']['UniProtACList'] = args.UniProtACList
        YAML_File['General informations']['CorrespondencesFile'] = args.CorrespondencesFile
        YAML_File['General informations']['OutputDirName'] = args.OutputDirName
        YAML_File['General informations']['ProjectDescription'] = args.ProjectDescription
        YAML_File['General informations']['IncludedPseudogenes'] = args.IncludedPseudogenes
        YAML_File = writtingYAML(
            YAML_File, outfile, 'General informations', ORDERBOX)

        # Clustering into Families
        step = ORDERBOX[2-gap]
        YAML_File[step] = {}
        YAML_File[step]['Coverage'] = args.Coverage
        YAML_File[step]['Identity'] = args.Identity
        addsAvancedSettings(
            YAML_File, step, args.MMseqsAdvancedSettings, common.getMMseqsDefaultSettings())
        YAML_File = writtingYAML(YAML_File, outfile, step, ORDERBOX)

        # Synteny Finder
        step = ORDERBOX[3-gap]
        YAML_File[step] = {}
        YAML_File[step]['WindowSize'] = args.WindowSize
        YAML_File[step]['SyntenyGap'] = args.SyntenyGap
        YAML_File[step]['SyntenyScoreCutoff'] = args.SyntenyScoreCutoff
        addsAvancedSettings(YAML_File, step, args.ClusteringAdvancedSettings,
                            common.getClusteringMethodsDefaultSettings())
        YAML_File = writtingYAML(YAML_File, outfile, step, ORDERBOX)

        # Data export
        step = ORDERBOX[4-gap]
        YAML_File[step] = {}
        YAML_File[step]['MetaDataFile'] = args.MetaDataFile
        YAML_File[step]['ClusteringMethod'] = args.ClusteringMethod
        YAML_File[step]['GroupingOnLabel'] = args.GroupingOnLabel
        YAML_File[step]['GroupingOnTaxonomy'] = args.GroupingOnTaxonomy
        YAML_File = writtingYAML(YAML_File, outfile, step, ORDERBOX)


def writtingYAML(YAML_to_write, outfile, block, ORDERBOX):
    '''
    YAML file writting.
    '''
    stream = yaml.dump(YAML_to_write, default_flow_style=False)
    if block == 'General informations':
        outfile.write('# ~~~~~ {}\n'.format('General informations'))
    elif block == ORDERBOX[1]:
        outfile.write('# ~~~~~ {}\n'.format(ORDERBOX[1]))
    else:
        outfile.write('# ~~~~~ Parameters for ' + block + '\n')
    outfile.write(stream)
    outfile.write('\n')
    return({})


def getBoxToResum(args):
    '''
    detection of settings modification
    make the association with box
    '''
    logger = logging.getLogger(getBoxToResum.__name__)
    boxToResume = None
    with open(common.global_dict['settingsFileName'], 'r') as file:
        oldSettings = yaml.load(file, Loader=yaml.FullLoader)
    newsettings = vars(args)
    advancedSettingsMMseqs = common.readYamlAdvancedSettingsFile(
        args.MMseqsAdvancedSettings, common.getMMseqsDefaultSettings())
    for settings in advancedSettingsMMseqs.values():
        for name, value in settings.items():
            newsettings[name] = value
    ClusteringAdvancedSettings = common.readYamlAdvancedSettingsFile(
        args.ClusteringAdvancedSettings, common.getClusteringMethodsDefaultSettings())
    for settings in ClusteringAdvancedSettings.values():
        for name, value in settings.items():
            newsettings[name] = value
    for boxKey in oldSettings:
        for parameterKey in oldSettings[boxKey]:
            if not newsettings[parameterKey] == oldSettings[boxKey][parameterKey]:
                if not boxToResume:
                    logger.info('New settings detected')
                    boxToResume = boxKey
                logger.info('New parameter for {}: {} -> {}'.format(parameterKey,
                                                                    oldSettings[boxKey][parameterKey], newsettings[parameterKey]))
    return boxToResume


def getLastBoxDone(ORDERBOX):
    '''
    Get name box executed.
    '''
    logger = logging.getLogger(getLastBoxDone.__name__)
    if os.path.isfile(common.global_dict['reportFileName']):
        with open(common.global_dict['reportFileName'], 'r') as file:
            nameBox = file.read()
        if not nameBox in ORDERBOX+['']:
            logger.critical('File of report wrong')
            exit(1)
        return nameBox
    return None


def inputFileModificationChecking(newInputFile, oldInputFile):
    '''
    Check that the input file has not been modified.
    '''
    logger = logging.getLogger(inputFileModificationChecking.__name__)
    error = common.checkFilledFile(newInputFile)
    if error:
        logger.error('The input file {} does not exist.'.format(newInputFile))
        exit(1)
    elif not filecmp.cmp(newInputFile, oldInputFile):
        logger.error(
            'The input file {} has been modified. Please start a new project.'.format(newInputFile))
        exit(1)


def inputFilesModificationChecking(inputI, inputII):
    '''
    Check that the input files has not been modified.
    '''
    if inputI:
        inputFileModificationChecking(
            inputI, common.global_dict['uniprotACListSaved'])
    if inputII:
        inputFileModificationChecking(
            inputII, common.global_dict['correspondingFileSaved'])


def resumptionFrom(nameBoxToResum, nameBoxDone, ORDERBOX):
    '''
    Returns the name of the box from which to resume the analysis.
    '''
    logger = logging.getLogger(resumptionFrom.__name__)
    if nameBoxDone:
        indexBoxDone = ORDERBOX.index(nameBoxDone)
    else:
        indexBoxDone = -1
    if nameBoxToResum:
        indexBoxToResum = ORDERBOX.index(nameBoxToResum)
    else:
        indexBoxToResum = -1
    if indexBoxDone == (len(ORDERBOX)-1):
        ''' Last analysis completed '''
        logger.debug('Last analysis completed')
        if indexBoxToResum == -1:
            ''' Not new analysis '''
            logger.debug('Not new analysis')
            return None
        else:
            ''' New Analysis '''
            logger.debug('New Analysis')
            return ORDERBOX[indexBoxToResum]
    else:
        ''' Last analysis not completed '''
        logger.debug('Last analysis not completed')
        if indexBoxToResum >= 0:
            ''' New Analysis '''
            logger.debug('New Analysis')
            if indexBoxToResum <= indexBoxDone+1:
                ''' New Analysis '''
                logger.debug('New Analysis')
                return ORDERBOX[indexBoxToResum]
            else:
                ''' Error Recovery '''
                logger.debug('Error Recovery')
                return ORDERBOX[indexBoxDone+1]
        else:
            ''' Error Recovery '''
            logger.debug('Error Recovery')
            return ORDERBOX[indexBoxDone+1]


def removeDownloadedINSDC():
    '''
    Remove dolowded INSDC files from GetINSDC step.
    '''
    boxName = common.global_dict['boxName']['GetINSDCFiles']
    dataDirectoryProcess = os.path.join(
        common.global_dict['dataDirectory'], boxName)
    filesExtension = common.global_dict['filesExtension']
    filesToRemoved = glob.glob(
        '{}/*.{}'.format(dataDirectoryProcess, filesExtension))
    for file in filesToRemoved:
        os.remove(file)


def boxesManager(runFrom, resultsDirectory, analysisNumber, ORDERBOX, args):
    '''
    Allows execution from a box.
    '''
    logger = logging.getLogger(boxesManager.__name__)
    boxesToRun = [ORDERBOX[i]
                  for i in range(ORDERBOX.index(runFrom), len(ORDERBOX))]
    logger.debug('Boxes to run: {}'.format(boxesToRun))
    for nameBox in boxesToRun:
        runBox(nameBox, resultsDirectory, analysisNumber, ORDERBOX, args)
        checkingAfterRun(nameBox)
        updateLastBoxDone(nameBox)


def updateLastBoxDone(nameBox):
    '''
    Updating report file.
    '''
    with open(common.global_dict['reportFileName'], 'w') as file:
        file.write(nameBox)


def getAnalysisNumber(outputDirName):
    '''
    Defines the number of analisys and defines the
    results directory name for the current analysis
    from the last completed analysis number.
    '''
    number = 0
    fileName = common.global_dict['lastAnalysisNumber']
    if os.path.isfile(fileName):
        with open(fileName, 'r') as file:
            number = int(file.read())
    analysisNumber = (number+1)
    resultsDirectory = '{}_Results_{}'.format(os.path.join(
        common.global_dict['workingDirectory'], outputDirName), analysisNumber)
    return (number+1), resultsDirectory


def checkingAfterRun(boxName):
    '''
    Checkpoints before executing box.
    '''
    logger = logging.getLogger(checkingAfterRun.__name__)
    pointsToCheck = common.global_dict['files']
    if boxName not in pointsToCheck.keys():
        logger.critical('No checkpoint for {}'.format(boxName))
        exit(1)
    if pointsToCheck[boxName]:
        error = False
        for fileName in pointsToCheck[boxName].keys():
            error = common.checkFilledFile(
                pointsToCheck[boxName][fileName], error)
        if error:
            logger.error(
                'Checkpoints not validated after {} execution'.format(boxName))
            exit(1)
    else:
        logger.error('Checkpoints is None for {}'.format(boxName))
        exit(1)


def mergeInputsII(outputStepI, inputII, inputMergedName, inputI):
    '''
    Merges the inputI and inputII.
    Copies of insdc files of inputII into the process directory of etINSDC step.
    '''
    logger = logging.getLogger(mergeInputsII.__name__)
    authorized_columns = common.definesAuthorizedColumns()
    mandatory_columns = common.definesMandatoryColumns()

    rowsI, _ = common.parseInputII(
        outputStepI, authorized_columns, mandatory_columns)
    rowsII, _ = common.parseInputII(
        inputII, authorized_columns, mandatory_columns)

    errors = False
    accessions = []
    for row in rowsI:
        accessions.append(row['protein_AC'])
        if 'taxon_ID' in rowsII[0].keys():
            row['taxon_ID'] = 'NA'
    for row in rowsII:
        if row['protein_AC'] in accessions:
            logger.error('{}: entry duplicated. A protein provided in {} is already related to this protein accession.\nPlease check redundancies between {} and {} files'.format(
                row['protein_AC'], inputI, outputStepI, inputII))
            errors = True
        else:
            # ----- pas nécessaire vu que parseInputII a déjà fait le travail sur les duplicats ----- ######
            accessions.append(row['protein_AC'])
            if not common.global_dict['inputIheader'] in row.keys():
                row[common.global_dict['inputIheader']] = 'NA'
            rowsI.append(row)
    if errors:
        logger.error('Merge of inputs failed')
        exit(1)
    else:
        with open(inputMergedName, 'w') as file:
            firstLine = True
            for row in rowsI:
                if firstLine:
                    headers = [header for header in row.keys()]
                    file.write('{}\n'.format('\t'.join(headers)))
                    firstLine = False
                line = ['' for i in row.values()]
                for header, value in row.items():
                    line[headers.index(header)] = value
                file.write('{}\n'.format('\t'.join(line)))
    return 0


def runBox(nameBox, resultsDirectory, analysisNumber, ORDERBOX, args):
    '''
    Running box.
    '''
    logger = logging.getLogger(runBox.__name__)
    if nameBox == 'GetINSDCFiles':
        netsyn_getINSDCFiles.run(args.UniProtACList)
    elif nameBox == 'ParseINSDCFiles_GetTaxonomy':
        logger = logging.getLogger('runParseINSDCFiles_GetTaxonomy')
        if args.UniProtACList and args.CorrespondencesFile:
            mergeInputsII(common.global_dict['files']['GetINSDCFiles']['inputClusteringStep'],
                          args.CorrespondencesFile,
                          common.global_dict['inputsMergedName'],
                          args.UniProtACList)
            inputBoxII = common.global_dict['inputsMergedName']
        elif args.CorrespondencesFile and not args.UniProtACList:
            inputBoxII = args.CorrespondencesFile
        else:
            inputBoxII = common.global_dict['files']['GetINSDCFiles']['inputClusteringStep']
        netsyn_parseINSDCFiles_GetTaxonomy.run(
            inputBoxII, args.IncludedPseudogenes)
        if 'GetINSDCFiles' in ORDERBOX:
            if not args.conserveDownloadedINSDC:
                logger.info('Removing downloded files in GetINSDC step')
                removeDownloadedINSDC()
            else:
                logger.info('Conserving downloded files in GetINSDC step')
    elif nameBox == 'ClusteringIntoFamilies':
        fasta = common.global_dict['files']['ParseINSDCFiles_GetTaxonomy']['faa']
        proteins = common.global_dict['files']['ParseINSDCFiles_GetTaxonomy']['proteins_1']
        netsyn_clusteringIntoFamilies.run(
            fasta, proteins, args.Identity, args.Coverage, args.MMseqsAdvancedSettings)
    elif nameBox == 'SyntenyFinder':
        logger = logging.getLogger('runSyntenyFinder')
        proteins = common.global_dict['files']['ClusteringIntoFamilies']['proteins_2']
        targets = common.global_dict['files']['ParseINSDCFiles_GetTaxonomy']['targets_2']
        netsyn_syntenyFinder.run(proteins, targets, args.WindowSize, args.SyntenyGap,
                                 args.SyntenyScoreCutoff, args.ClusteringAdvancedSettings)
    elif nameBox == 'DataExport':
        nodesFile = common.global_dict['files']['SyntenyFinder']['nodes']
        edgesFile = common.global_dict['files']['SyntenyFinder']['edges']
        metricsFile = common.global_dict['files']['SyntenyFinder']['metrics']
        organimsFile = common.global_dict['files']['ParseINSDCFiles_GetTaxonomy']['organisms_2']
        proteinsFile = common.global_dict['files']['SyntenyFinder']['proteins']
        netsyn_dataExport.run(nodesFile,
                              edgesFile,
                              organimsFile,
                              proteinsFile,
                              metricsFile,
                              args.MetaDataFile,
                              args.GroupingOnLabel,
                              args.GroupingOnTaxonomy,
                              args.ClusteringMethod)
    elif nameBox == 'EndNetSynAnalysis':
        EndNetSynAnalysis.run(resultsDirectory, analysisNumber, ORDERBOX)
    else:
        logger.critical('No execution planned for {}'.format(nameBox))
        exit(1)


def main():
    ######################
    # Parse command line #
    ######################
    args = argumentsParser()
    ##########
    # Logger #
    ##########
    common.parametersLogger(args)
    logger = logging.getLogger('main')
    #########################
    # Dependancies checking #
    #########################
    common.dependenciesChecking()
    ##############
    # Constantes #
    ##############
    args.OutputDirName = args.OutputDirName.strip('/')
    common.constantsInitialization(
        args.OutputDirName, args.UniProtACList, args.CorrespondencesFile)
    ORDERBOX = orderDefinition(args.UniProtACList, args.CorrespondencesFile)
    logger.debug('ORDERBOX: {}'.format(ORDERBOX))
    ############
    # Switches #
    ############
    if args.newProject and os.path.isdir(common.global_dict['workingDirectory']):
        shutil.rmtree(common.global_dict['workingDirectory'])
    checkSameVersion(common.global_dict['versionFileName'])
    if not os.path.isfile(common.global_dict['settingsFileName']):
        if os.path.isdir(common.global_dict['workingDirectory']):
            logger.error('{} already existing directory and not dedicated to NetSyn project. Please change your project name or use --newProject option'.format(
                common.global_dict['workingDirectory']))
            logger.error(
                'Project name already used. Please change your project name or use --newProject.')
            exit(1)
        '''Pas de dernier setting == premiere analyse'''
        logger.info('New NetSyn project: {}'.format(
            common.global_dict['workingDirectory']))
        os.mkdir(common.global_dict['workingDirectory'])
        os.mkdir(common.global_dict['dataDirectory'])
        analysisNumber, resultsDirectory = getAnalysisNumber(
            args.OutputDirName)
        common.filesNameInitialization(
            resultsDirectory, args.OutputDirName, analysisNumber)
        try:
            saveInputFiles(args.UniProtACList, args.CorrespondencesFile)
        except:
            shutil.rmtree(common.global_dict['workingDirectory'])
            exit(1)
        runFromBox = ORDERBOX[0]
        nameBoxToResum = None
        nameBoxDone = None
    else:
        '''Fichier de setting == analyse deja lance'''
        nameBoxToResum = getBoxToResum(args)
        if nameBoxToResum == 'General informations':
            logger.error(
                'Already existing project. Please start a new project')
            exit(1)
        else:
            inputFilesModificationChecking(
                args.UniProtACList, args.CorrespondencesFile)
            logger.info('NetSyn project recovery: {}'.format(
                common.global_dict['workingDirectory']))
        nameBoxDone = getLastBoxDone(ORDERBOX)
        runFromBox = resumptionFrom(nameBoxToResum, nameBoxDone, ORDERBOX)
        analysisNumber, resultsDirectory = getAnalysisNumber(
            args.OutputDirName)
        common.filesNameInitialization(
            resultsDirectory, args.OutputDirName, analysisNumber)
    createYamlSettingsFile(args, ORDERBOX)
    if runFromBox:
        '''Execution a partir de la boite x'''
        runFromBoxIndex = ORDERBOX.index(runFromBox)
        updateLastBoxDone(ORDERBOX[runFromBoxIndex-1]
                          if runFromBoxIndex-1 >= 0 else '')
        logger.debug('Box from which resume: {}'.format(nameBoxToResum))
        logger.debug('Last box done: {}'.format(nameBoxDone))
        logger.debug('Run from: {}'.format(runFromBox))
        writeCurrentVersion(
            common.global_dict['version'], common.global_dict['versionFileName'])
        boxesManager(runFromBox, resultsDirectory,
                     analysisNumber, ORDERBOX, args)
    else:
        '''Nous avons deja des resultats'''
        logger.debug('Box from which resume: {}'.format(nameBoxToResum))
        logger.debug('Last box done: {}'.format(nameBoxDone))
        logger.debug('Run from: {}'.format(runFromBox))
        logger.info('This analysis already exists')


####################
# Script execution #
####################
if __name__ == '__main__':
    main()
