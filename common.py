##########
# Import #
##########
import __main__ as namespace
import json
import pickle
import os
import logging
import urllib3
import subprocess
import errno

#############
# Functions #
#############
def unPetitBonjourPourredonnerLeMoral(msg=None):
    '''
    Fonction de desespoire pour les developpeur en quete de debogage...
    '''
    print("=D I'M HAPPY!!!")
    if msg:
        print(msg)

def widowsSizePossibilities(minSize, maxSize):
    return range(minSize, maxSize+2, 2)

def dependanciesChecking():
    logger = logging.getLogger('{}.{}'.format(checkFilledFile.__module__, checkFilledFile.__name__))
    try:
        devnull = open(os.devnull)
        subprocess.Popen('mmseqs', stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            logger.error('mmseqs not found. Please check its installation.')
            exit(1)

def checkFilledFile(fileName, error=False):
    '''
    Checking if file existing and not empty.
    '''
    logger = logging.getLogger('{}.{}'.format(checkFilledFile.__module__, checkFilledFile.__name__))
    if not os.path.isfile(fileName):
        error = True
        logger.error('{} missing.'.format(fileName))
    elif os.path.getsize(fileName) == 0:
        error = True
        logger.error('{} empty.'.format(fileName))
    return error

def httpRequest(poolManager,method, url):
    '''
    Return http request result.
    '''
    logger = logging.getLogger('{}.{}'.format(httpRequest.__module__, httpRequest.__name__))
    retry = urllib3.util.Retry(read=5, backoff_factor=2)
    try:
        res = poolManager.request(method, url, retries=retry)
    except urllib3.exceptions.NewConnectionError:
        logger.error('Connection failed.')
        exit(1)
    return res

def constantsInitialization(outputDirName, inputFile):
    '''
    Initialization of constants.
    Calling from netsyn or BoxManager modules.
    '''
    global_dict['workingDirectory'] = outputDirName
    global_dict['dataDirectory'] = '{}/data'.format(outputDirName) #TMPDIRECTORY
    global_dict['settingsFileName'] = '{}/{}'.format(outputDirName, '.lastSettings.yml') #SETTINGSFILENAME
    global_dict['reportFileName'] = '{}/{}'.format(outputDirName, '.report') #REPORTFILENAME
    global_dict['versionFileName'] = '{}/{}'.format(outputDirName, '.version') #REPORTFILENAME
    global_dict['inputFileSaved'] = '{}/{}'.format(outputDirName, os.path.basename(inputFile)) #INPUTLIST

def filesNameInitialization(resultsDirectory, outputDirName, analysisNumber):
    global_dict['files'] = {
        global_dict['boxName']['GetINSDCFiles'] : {
            'inputClusteringStep' : '{}/{}/inputClusteringIntoFamiliesStep.tsv'.format(global_dict['dataDirectory'], global_dict['boxName']['GetINSDCFiles'])
        },
        global_dict['boxName']['ClusteringIntoFamilies'] : {
            'faa' : '{}/{}/MMseqs2_run.faa'.format(global_dict['dataDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'contigs' : '{}/{}/contigs.pickle'.format(global_dict['dataDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'genomicContexts' : '{}/{}/genomicContexts.pickle'.format(global_dict['dataDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'lineage' : '{}/{}/taxonomicLineage.pickle'.format(global_dict['dataDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'targets' : '{}/{}/targets_list.pickle'.format(global_dict['dataDirectory'], global_dict['boxName']['ClusteringIntoFamilies'])
        },
        global_dict['boxName']['SyntenyFinder'] : {
            'nodes': '{}/{}/nodes_list.pickle'.format(global_dict['dataDirectory'], global_dict['boxName']['SyntenyFinder']),
            'edges': '{}/{}/edges_list.pickle'.format(global_dict['dataDirectory'], global_dict['boxName']['SyntenyFinder'])
        },
        global_dict['boxName']['DataExport'] : {
            'graphML' : '{}/{}_Results_{}.graphML'.format(resultsDirectory, outputDirName, analysisNumber),
            'nodes': '{}/{}_Results_{}.nodes.json'.format(resultsDirectory, outputDirName, analysisNumber),
            'edges': '{}/{}_Results_{}.edges.json'.format(resultsDirectory, outputDirName, analysisNumber),
            #'html' : '{}/{}_Results_{}.html'.format(resultsDirectory, outputDirName, analysisNumber),
            #'settings' : '{}/{}_Settings_{}.yaml'.format(resultsDirectory, outputDirName, analysisNumber)
        }
    }

def write_pickle(dictionary, output):
    '''
    '''
    #logger = logging.getLogger('{}.{}'.format(write_pickle.__module__, write_pickle.__name__))
    with open(output, 'wb') as pickleFile:
        pickle.dump(dictionary, pickleFile)
    return 0

def write_json(dictionary, output):
    '''
    '''
    #logger = logging.getLogger('{}.{}'.format(write_json.__module__, write_json.__name__))
    with open(output, 'w') as jsonFile:
        json.dump(dictionary, jsonFile, indent=4)
    return 0

def read_pickle(input):
    '''
    '''
    with open(input, 'rb') as file:
        return pickle.load(file)

def read_file(input):
    '''
    '''
    with open(input, 'r') as file:
        return file.readlines()

def parametersLogger(args):
    '''
    Logger setting.
    '''
    if not args.log_file:
        cyanColor = '\033[36m'
        yellowColor = '\033[33m'
        endColor = '\033[0m'
    else:
        cyanColor = ''
        yellowColor = ''
        endColor = ''
    logging_std_format = '{}[%(levelname)s]{} %(message)s'.format(yellowColor, endColor)
    logging_debug_format = '{}%(asctime)s {}[%(levelname)s]{} [%(threadName)s - %(name)s]{} %(message)s'.format(cyanColor, yellowColor, cyanColor, endColor)
    log_level = args.log_level.upper()
    if (log_level == 'DEBUG'):
        logging_std_format = logging_debug_format
    logging_datefmt = '%Y/%m/%d - %H:%M:%S'
    if (args.log_file != None):
        logging.basicConfig(format = logging_std_format,
                             datefmt = logging_datefmt,
                             filename = args.log_file,
                             filemode = 'w',
                             level = log_level)
    else:
        logging.basicConfig(format = logging_std_format,
                             datefmt = logging_datefmt,
                             level = log_level)
#########################
# Constantes definition #
#########################
global_dict = {
    'version': '0.0.1',
    'defaultValue': 'NA',
    'maxGCSize': 11, #MAXGCSIZE
    'minGCSize': 3,
    'boxName' : {
        'GetINSDCFiles' : 'GetINSDCFiles',
        'ClusteringIntoFamilies' :'ClusteringIntoFamilies',
        'SyntenyFinder' : 'SyntenyFinder',
        'DataExport' : 'DataExport'
    },
    'inputIIheaders' : [
        'UniProt_AC',
        'protein_AC',
        'protein_AC_field',
        'nucleic_AC',
        'nucleic_File_Format',
        'nucleic_File_Name'
    ],
    'desired_ranks_lneage' : {
        'superkingdom' : 1,
        'phylum' : 5,
        'class' : 8,
        'order' : 13,
        'family' : 17,
        'genus' : 21,
        'species' : 26
    }
}


##############################
# Add variables to namespace #
##############################
for variable,value in global_dict.items():
    setattr(namespace, variable, value)
