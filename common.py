##########
# Import #
##########
import __main__ as namespace
import json
import pickle
import os
import logging
import urllib3

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

def constantsInitialization(projectName, inputFile):
    '''
    Initialization of constants.
    Calling from netsyn or BoxManager modules.
    '''
    global_dict['workingDirectory'] = projectName
    global_dict['tmpDirectory'] = '{}/TMP'.format(projectName) #TMPDIRECTORY
    global_dict['settingsFileName'] = '{}/{}'.format(projectName, '.lastSettings.yml') #SETTINGSFILENAME
    global_dict['reportFileName'] = '{}/{}'.format(projectName, '.report.yml') #REPORTFILENAME
    global_dict['inputFileSaved'] = '{}/{}'.format(projectName, inputFile) #INPUTLIST

def filesNameInitialization(resultsDirectory, projectName, analysisNumber):
    global_dict['files'] = {
        global_dict['boxName']['GetINSDCFiles'] : {
            'inputClusteringStep' : '{}/{}/inputClusteringIntoFamiliesStep.tsv'.format(global_dict['tmpDirectory'], global_dict['boxName']['GetINSDCFiles'])
        },
        global_dict['boxName']['ClusteringIntoFamilies'] : {
            'faa' : '{}/{}/MMseqs2_run.faa'.format(global_dict['tmpDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'contigs' : '{}/{}/contigs.pickle'.format(global_dict['tmpDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'genomicContexts' : '{}/{}/genomicContexts.pickle'.format(global_dict['tmpDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'lineage' : '{}/{}/taxonomicLineage.pickle'.format(global_dict['tmpDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'targets' : '{}/{}/targets_list.pickle'.format(global_dict['tmpDirectory'], global_dict['boxName']['ClusteringIntoFamilies'])
        },
        global_dict['boxName']['SyntenyFinder'] : {
            'nodes': '{}/{}/nodes_list.pickle'.format(global_dict['tmpDirectory'], global_dict['boxName']['SyntenyFinder']),
            'edges': '{}/{}/edges_list.pickle'.format(global_dict['tmpDirectory'], global_dict['boxName']['SyntenyFinder'])
        },
        global_dict['boxName']['DataExport'] : {
            'graphML' : '{}/{}_Results_{}.graphML'.format(resultsDirectory, projectName, analysisNumber),
            'nodes': '{}/{}_Results_{}.nodes.json'.format(resultsDirectory, projectName, analysisNumber),
            'edges': '{}/{}_Results_{}.edges.json'.format(resultsDirectory, projectName, analysisNumber),
            #'html' : '{}/{}_Results_{}.html'.format(resultsDirectory, projectName, analysisNumber),
            #'settings' : '{}/{}_Settings_{}.yaml'.format(resultsDirectory, projectName, analysisNumber)
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

#########################
# Constantes definition #
#########################
global_dict = {
    'defaultValue': 'NA',
    'maxGCSize': 11, #MAXGCSIZE
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
