##########
# Import #
##########
import __main__ as namespace
import json
import pickle
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


def constantsInitialiszation(projectName, inputFile,):
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
        return file.read()

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
    ]
}


##############################
# Add variables to namespace #
##############################
for variable,value in global_dict.items():
    setattr(namespace, variable, value)
