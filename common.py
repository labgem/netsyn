##########
# Import #
##########
import __main__ as namespace
import json
import os
import logging
import urllib3
import subprocess
import errno
import re

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

def parseInputI(filename): # Fonction a deplacer dans tools ??
    '''
    Input file parssing.
    '''
    logger = logging.getLogger('{}.{}'.format(parseInputI.__module__, parseInputI.__name__))
    with open(filename, 'r') as file:
        error = False
        firstLine = True
        accessions = []
        seps = [' ','\t',',',';']
        for line in file:
            for sep in seps:
                if len(line.split(sep)) > 1:
                    logger.error('Input invalidated: unauthorized character {}'.format(seps))
                    exit(1)
            value = line.rstrip()
            if firstLine:
                header = value
                firstLine = False
            elif not value in accessions:
                accessions.append(value)
            else:
                error = True
                logger.error('{}: Entry duplicated'.format(value))
    if error:
        logger.error('Input invalidated.')
        exit(1)
    return header, accessions

def checkInputHeaders(errors, mandatory_columns, headers):
    ''' check if all mandatory columns are provided
    '''
    logger = logging.getLogger('{}.{}'.format(checkInputHeaders.__module__, checkInputHeaders.__name__))
    for mandatory_column in mandatory_columns:
        if not mandatory_column in headers.values():
            logger.error('{}: Missing column'.format(mandatory_column))
            errors = True
    return errors

def parseInputII(fname, authorized_columns, mandatory_columns):
    ''' create a list of dictionaries
        get from input file (fname) information by line
        every line is a dictionary stored in a list
    '''
    logger = logging.getLogger('{}.{}'.format(parseInputII.__module__, parseInputII.__name__))
    first_line = True
    errors = False
    rows = []
    line_number = 0
    with open(fname, 'r') as file:
        accessions = []
        for line in file:
            if first_line:
                headers = {}
                for index, header in enumerate(line.split('\t')):
                    header = header.replace('\r\n', '').replace('\n', '') # header.strip() ???
                    p = re.compile(r'(?:{})'.format('|'.join(authorized_columns)))
                    if not p.search(header):
                        logger.error('{}: Column name not valid'.format(header))
                        errors = True
                    if header in headers.values():
                        logger.error('{}: Duplicated column'.format(header))
                        errors = True
                    headers[index] = header
                errors = checkInputHeaders(errors, mandatory_columns, headers)
                first_line = False
            else:
                line_number += 1
                row = {}
                for index, column in enumerate(line.strip().split('\t')):
                    if column == '':
                        logger.error('Empty field: line {}, column "{}"'.format(line_number, headers[index]))
                        errors = True
                    elif headers[index] == proteinACHeader:
                        if column in accessions:
                            logger.error('{}: duplicated entry'.format(column))
                            errors = True
                        else:
                            accessions.append(column)
                    elif headers[index] == global_dict['inputIIheaders'][global_dict['inputIIheaders'].index('nucleic_File_Path')]:
                        errors = checkFilledFile(column, errors)

                    row[headers[index]] = column.replace('\r\n', '').replace('\n', '') # header.strip()
                rows.append(row)
    if errors:
        logger.error('Invalid input: {}'.format(fname))
        exit(1)
    return rows, list(headers.values())

def definesAuthorizedColumns():
    '''
    Defines the authorized columns.
    '''
    return global_dict['inputIIheaders'] + ['taxon_ID']

def definesMandatoryColumns():
    '''
    Defines the mandatory columns.
    '''
    mandatory_columns = list(global_dict['inputIIheaders'])
    mandatory_columns.remove('UniProt_AC')
    return mandatory_columns

def widowsSizePossibilities(minSize, maxSize):
    return range(minSize, maxSize+2, 2)

def dependenciesChecking():
    logger = logging.getLogger('{}.{}'.format(dependenciesChecking.__module__, dependenciesChecking.__name__))
    try:
        devnull = open(os.devnull)
        subprocess.Popen('mmseqs', stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            logger.error('mmseqs not found. Please check its installation')
            exit(1)

def checkFilledFile(fileName, error=False):
    '''
    Checking if file existing and not empty.
    '''
    logger = logging.getLogger('{}.{}'.format(checkFilledFile.__module__, checkFilledFile.__name__))
    if not os.path.isfile(fileName):
        error = True
        logger.error('{} missing'.format(fileName))
    elif os.path.getsize(fileName) == 0:
        error = True
        logger.error('{} empty'.format(fileName))
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

def constantsInitialization(outputDirName, uniprotACList, correspondingFile):
    '''
    Initialization of constants.
    Calling from netsyn or BoxManager modules.
    '''
    global_dict['workingDirectory'] = outputDirName
    global_dict['dataDirectory'] = os.path.join(outputDirName, 'data')
    global_dict['inputsMergedName'] = os.path.join(outputDirName, 'inputsMerged.tsv')
    global_dict['settingsFileName'] = os.path.join(outputDirName, '.lastSettings.yml')
    global_dict['reportFileName'] = os.path.join(outputDirName, '.report')
    global_dict['versionFileName'] = os.path.join(outputDirName, '.version')
    global_dict['lastAnalysisNumber'] = os.path.join(outputDirName, '.lastAnalysisNumber')
    if uniprotACList:
        global_dict['uniprotACListSaved'] = os.path.join(outputDirName, os.path.basename(uniprotACList))
    if correspondingFile:
        global_dict['correspondingFileSaved'] = os.path.join(outputDirName, os.path.basename(correspondingFile))

def filesNameInitialization(resultsDirectory, outputDirName, analysisNumber):
    global_dict['files'] = {
        global_dict['boxName']['GetINSDCFiles'] : {
            'inputClusteringStep' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['GetINSDCFiles'], 'inputClusteringIntoFamiliesStep.tsv')
        },
        global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'] : {
            'proteins_1' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'proteins_parsingStep.json'),
            'organisms_1' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'organisms_parsingStep.json'),
            'organisms_2' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'organisms_taxonomyStep.json'),
            'targets_1' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'targets_parsingStep.json'),
            'targets_2' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'targets_taxonomyStep.json'),
            'faa' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'MMseqs2_run.faa')
        },
        global_dict['boxName']['ClusteringIntoFamilies'] : {
            'proteins_2' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ClusteringIntoFamilies'], 'proteins_familiesStep.json')
        },
        global_dict['boxName']['SyntenyFinder'] : {
            'nodes': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['SyntenyFinder'], 'nodes_list.json'),
            'edges': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['SyntenyFinder'], 'edges_list.json'),
            'proteins': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['SyntenyFinder'], 'proteins_syntenyStep.json')
        },
        global_dict['boxName']['DataExport'] : {
            'graphML' : '{}_Results_{}.graphML'.format(os.path.join(resultsDirectory, outputDirName), analysisNumber),
            'html' : '{}_Results_{}.html'.format(os.path.join(resultsDirectory, outputDirName), analysisNumber),
        },
        global_dict['boxName']['EndNetSynAnalysis'] : {
            'settings' : '{}_Settings_{}.yaml'.format(os.path.join(resultsDirectory, outputDirName), analysisNumber)
        }
    }

def write_json(dictionary, output):
    '''
    '''
    with open(output, 'w') as jsonFile:
        json.dump(dictionary, jsonFile, indent=4)
    return 0

def readJSON(nameFile):
    '''
    '''
    with open(nameFile, 'r') as file:
        return json.load(file)

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
inputIheader = 'UniProt_AC'
proteinACHeader = 'protein_AC'
global_dict = {
    'version': '0.0.3',
    'defaultValue': 'NA',
    'maxGCSize': 11, #MAXGCSIZE
    'minGCSize': 3,
    'filesExtension': 'embl',
    'boxName': {
        'GetINSDCFiles': 'GetINSDCFiles',
        'ParseINSDCFiles_GetTaxonomy': 'ParseINSDCFiles_GetTaxonomy',
        'ClusteringIntoFamilies': 'ClusteringIntoFamilies',
        'SyntenyFinder': 'SyntenyFinder',
        'DataExport': 'DataExport',
        'EndNetSynAnalysis' : 'EndNetSynAnalysis'
    },
    'inputIheader': inputIheader,
    'proteinACHeader': proteinACHeader,
    'inputIIheaders': [
        inputIheader,
        proteinACHeader,
        'protein_AC_field',
        'nucleic_AC',
        'nucleic_File_Format',
        'nucleic_File_Path'
    ],
    'desired_ranks_lneage': {
        'superkingdom': 1,
        'phylum': 5,
        'class': 8,
        'order': 13,
        'family': 17,
        'genus': 21,
        'species': 26
    },
    'metadataMadatoryColumn': [
        'accession_type',
        'accession'
    ],
    'metadataAccessionAuthorized': [
        inputIheader,
        proteinACHeader
    ]
}

##############################
# Add variables to namespace #
##############################
for variable,value in global_dict.items():
    setattr(namespace, variable, value)
