##########
# Import #
##########
import __main__ as namespace
from netsyn import __version__

import json
import jsonschema
import os
import logging
import urllib3
import subprocess
import errno
import re
import yaml

#############
# Functions #
#############
def getMMseqsDefaultSettings():
    '''
    Defines MMseqs2 default Seetings.
    '''
    return {
        'MMseqs advanced settings': {
            'MMseqs_cov-mode': 1,
            'MMseqs_cluster-mode': 1,
            'MMseqs_kmer-per-seq': 80,
            'MMseqs_max-seqs': 300,
            'MMseqs_single-step-clustering': 'false',
            'MMseqs_threads': 4
        }
    }

def getClusteringMethodsDefaultSettings():
    '''
    Defines clustering methods default Seetings.
    '''
    return {
        global_dict['MCL']: {
            'MCL_inflation': 2,
            'MCL_expansion': 2,
            'MCL_iterations': 1000
        },
        global_dict['WalkTrap']: {
            'walktrap_step': 4
        },
        global_dict['Infomap']: {
            'infomap_trials': 10
        }
    }

def readYamlAdvancedSettingsFile(yamlFileName, defaultSettings):
    '''
    Read the advanced settings file and compare to default settings.
    '''
    logger = logging.getLogger('{}.{}'.format(readYamlAdvancedSettingsFile.__module__, readYamlAdvancedSettingsFile.__name__))
    error = False
    if yamlFileName:
        with open(yamlFileName, 'r') as file:
            content = yaml.load(file, Loader=yaml.BaseLoader)
        for name, settings in content.items():
            for setting, value in settings.items():
                if name not in defaultSettings.keys():
                    logger.error('{} is a settings name not allowed.'.format(name))
                    logger.error('Name allowed: {}'.format(', '.join(defaultSettings.keys())))
                    error = True
                elif setting not in defaultSettings[name].keys():
                    logger.error('{} is not allowed as mmseqs setting.'.format(setting))
                    logger.error('Settings allowed: {}'.format(', '.join(defaultSettings[name].keys())))
                    error = True
                else:
                    try:
                        defaultSettings[name][setting] = int(value)
                    except:
                        defaultSettings[name][setting] = value

    if error:
        logger.error('The {} file is improper.'.format(yamlFileName))
        exit(1)
    return defaultSettings


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
    if checkFilledFile(fname):
        logger.error('Please make sure that {} file is in the appropriate repertory'.format(fname))
        exit(1)
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
    retry = urllib3.util.Retry(read=5, connect=5, backoff_factor=0.5, status_forcelist=set([504]))
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
    global_dict['dataDirectory'] = os.path.join(outputDirName, '.netsyn')
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
            'inputClusteringStep' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['GetINSDCFiles'], 'inputClusteringIntoFamiliesStep.tsv'),
            'report': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['GetINSDCFiles'], 'report.txt')
        },
        global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'] : {
            'proteins_1' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'proteins_parsingStep.json'),
            'organisms_1' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'organisms_parsingStep.json'),
            'organisms_2' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'organisms_taxonomyStep.json'),
            'targets_1' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'targets_parsingStep.json'),
            'targets_2' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'targets_taxonomyStep.json'),
            'faa' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'multifasta.faa'),
            'report': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ParseINSDCFiles_GetTaxonomy'], 'report.txt')
        },
        global_dict['boxName']['ClusteringIntoFamilies'] : {
            'proteins_2' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ClusteringIntoFamilies'], 'proteins_familiesStep.json'),
            'families' : os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ClusteringIntoFamilies'], 'families.tsv'),
            'report': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['ClusteringIntoFamilies'], 'report.txt')
        },
        global_dict['boxName']['SyntenyFinder'] : {
            'nodes': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['SyntenyFinder'], 'nodes_list.json'),
            'edges': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['SyntenyFinder'], 'edges_list.json'),
            'proteins': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['SyntenyFinder'], 'proteins_syntenyStep.json'),
            'report': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['SyntenyFinder'], 'report.txt')
        },
        global_dict['boxName']['DataExport'] : {
            'graphML': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['DataExport'], 'Results.graphML'),
            'html': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['DataExport'], 'Results.html'),
            'report': os.path.join(global_dict['dataDirectory'], global_dict['boxName']['DataExport'], 'report.txt')
        },
        global_dict['boxName']['EndNetSynAnalysis'] : {
            'graphML': '{}_Results_{}.graphML'.format(os.path.join(resultsDirectory, outputDirName), analysisNumber),
            'html': '{}_Results_{}.html'.format(os.path.join(resultsDirectory, outputDirName), analysisNumber),
            'report': '{}_Report_{}.txt'.format(os.path.join(resultsDirectory, outputDirName), analysisNumber),
            'settings' : '{}_Settings_{}.yaml'.format(os.path.join(resultsDirectory, outputDirName), analysisNumber)
        }
    }
    global_dict['synthesisDataExport'] = os.path.join(global_dict['dataDirectory'], global_dict['boxName']['DataExport'], 'ClusteringsSynthesis')
    global_dict['synthesisEndNetSynAnalysis'] = '{}_ClusteringsSynthesis_{}'.format(os.path.join(resultsDirectory, outputDirName), analysisNumber)

def write_json(dictionary, output):
    '''
    '''
    with open(output, 'w') as jsonFile:
        json.dump(dictionary, jsonFile)#, indent=4)
    return 0

def getEdgesListStepschema():
    return {
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "source": {"type": "number"},
                "target": {"type": "number"},
                "proteins_idx_source": {
                    "type": "array",
                    "items": {"type": "string"},
                },
                "proteins_idx_target": {
                    "type": "array",
                    "items": {"type": "string"},
                },
                "weight": {"type": "number"},
            }, "required": [
                    "source",
                    "target",
                    "proteins_idx_source",
                    "proteins_idx_target",
                    "weight"
             ],

        },
    }

def getNodesListStepschema():
    return {
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "protein_idx": {"type": "number"},
                "id": {"type": "string"},
                "UniProt_AC": {"type": "string"},
                "protein_AC": {"type": "string"},
                "context": {
                    "type": "array",
                    "items": {"type": "string"},
                },
                "context_idx": {
                    "type": "array",
                    "items": {"type": "string"},
                },
                "organism_id": {"type": "number"},
                "organism_idx": {"type": "number"},
                "clusterings": {
                    "type": "object",
                },
                "families": {
                    "type": "array",
                    "items": {"type": "number"},
                },
                "Size": {"type": "number"},
            },
            "required": [
                "protein_idx",
                "id",
                "UniProt_AC",
                "protein_AC",
                "context",
                "context_idx",
                "organism_id",
                "organism_idx",
                "clusterings",
                "families",
                "Size",
            ],

        },
    }

def getTargetsTaxonomyStepschema():
    return {
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "id": {"type": "string"},
                "protein_idx": {"type": "string"},
                "organism_id": {"type": "number"},
                "context": {
                    "type": "array",
                    "items": {"type": "string"},
                },
                "context_idx": {
                    "type": "array",
                    "items": {"type": "string"},
                },
                "UniProt_AC": {"type": "string"},
                "protein_AC": {"type": "string"},
                "organism_idx": {"type": "number"},
            }, "required" : [
                    "id",
                    "protein_idx",
                    "organism_id",
                    "context",
                    "context_idx",
                    "UniProt_AC",
                    "protein_AC",
                    "organism_idx"
            ],
        },
    }

def getOrganismsTaxonomyStepschema():
    return {
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "id": {"type" : "number"},
                "name": {"type" : "string"},
                "strain": {"type" : "string"},
                "taxon_id": {"type" : "string"},
                "lineage": {
                    "type" : "array",
                    "items": {
                        "type" : "object",
                        "properties": {
                            "rank": {"type" : "string"},
                            "scienticName": {"type" : "string"},
                            "tax_id": {"type" : "string"},
                            "level": {"type" : "number"},
                        },
                        "required": [
                            "rank",
                            "scientificName",
                            "tax_id",
                            "level"
                        ]
                    },
                },
            },
            "required": [
                "id",
                "name",
                "strain",
                "taxon_id",
                "lineage"
            ]
        }
    }

def getProteinsParsingStepSchema():
    return {
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "id": {"type" : "string"},
                "protein_AC": {"type" : "string"},
                "begin": {"type" : "number"},
                "end": {"type" : "number"},
                "strand": {"type" : "string"},
                "products": {"type" : "string"},
                "ec_numbers": {"type" : "string"},
                "UniProt_AC": {"type" : "string"},
                "gene_names": {"type" : "string"},
                "locus_tag": {"type" : "string"},
                "targets": {"type" : "array", "items": {"type": "string"}},
                "targets_idx": {"type" : "array", "items": {"type": "string"}}
            }, "required" : [
                    "id",
                    "protein_AC",
                    "begin",
                    "end",
                    "strand",
                    "products",
                    "ec_numbers",
                    "UniProt_AC",
                    "gene_names",
                    "locus_tag",
                    "targets",
                    "targets_idx"
                ]
        },
    }

def getProteinsFamiliesStepSchema():
    schema = getProteinsParsingStepSchema()
    schema["items"]["properties"]["family"] = {"type" : "number"}
    schema["items"]["required"].append("family")
    return schema

def getProteinsSyntenyStepSchema():
    return getProteinsFamiliesStepSchema()

def readJSON(nameFile, schema):
    '''
    '''
    logger = logging.getLogger('{}.{}'.format(readJSON.__module__, readJSON.__name__))
    with open(nameFile, 'r') as file:
        data = json.load(file)
    # if type(schema) == type([]):
    #     schema = schema[0]
    #     if not type(data) == type([]):
    #         logger.error('Data from {} file must be into a list.'.format(nameFile))
    #         exit(0)
    #     for d in data:
    #         validate(d, schema)
    # else:
    #     validate(data, schema)
    if schema:
        try:
            jsonschema.validate(data, schema)
        except 	jsonschema.exceptions.UnknownType    as e:
            logger.error(e)
            exit(1)

    return data

def read_file(input):
    '''
    '''
    with open(input, 'r') as file:
        return file.readlines()

def parametersLogger(args):
    '''
    Logger setting.
    '''
    if not args.logFile:
        cyanColor = '\033[36m'
        yellowColor = '\033[33m'
        endColor = '\033[0m'
    else:
        cyanColor = ''
        yellowColor = ''
        endColor = ''
    logging_std_format = '{}[%(levelname)s]{} %(message)s'.format(yellowColor, endColor)
    logging_debug_format = '{}%(asctime)s {}[%(levelname)s]{} [%(threadName)s - %(name)s]{} %(message)s'.format(cyanColor, yellowColor, cyanColor, endColor)
    logLevel = args.logLevel.upper()
    if (logLevel == 'DEBUG'):
        logging_std_format = logging_debug_format
    logging_datefmt = '%Y/%m/%d - %H:%M:%S'
    if (args.logFile != None):
        logging.basicConfig(format = logging_std_format,
                             datefmt = logging_datefmt,
                             filename = args.logFile,
                             filemode = 'w',
                             level = logLevel)
    else:
        logging.basicConfig(format = logging_std_format,
                             datefmt = logging_datefmt,
                             level = logLevel)

def reportingFormat(logger, boxName, messages):
    '''
    Formats the reporting.
    '''
    if (messages):
        msg = '# {} reporting #'.format(boxName)
        contentFile = []
        logger.info('#'*len(msg))
        contentFile.append('#'*len(msg))
        logger.info(msg)
        contentFile.append(msg)
        logger.info('#'*len(msg))
        contentFile.append('#'*len(msg))
        for msg in messages:
            logger.info(msg)
            contentFile.append(msg)
        logger.info('#'*len(msg))
        contentFile.append('\n')
        with open(global_dict['files'][boxName]['report'], 'w') as file:
            file.write('\n'.join(contentFile))
    else:
        with open(global_dict['files'][boxName]['report'], 'w') as file:
            file.write('*')

#########################
# Constantes definition #
#########################
inputIheader = 'UniProt_AC'
proteinACHeader = 'protein_AC'
resourcesDir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'resources')
global_dict = {
    'version': __version__,
    'defaultValue': 'NA',
    'maxGCSize': 11, #MAXGCSIZE
    'minGCSize': 3,
    'sscDefault': 3.0,
    'filesExtension': 'embl',
    'formatOfFilesToParse': ['embl', 'gb', 'genbank'],
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
    'desired_ranks_lineage': {
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
    ],
    'graphML_edgesAttributes': [
        'protein_AC',
        'products',
        'ec_numbers',
        'UniProt_AC',
        'gene_names',
        'locus_tag'
    ],
    'Clustering_Methods': [
        'Infomap',
        'Louvain',
        'MCL',
        'WalkTrap'
    ],
    'MCL': 'MCL advanced settings',
    'WalkTrap': 'WalkTrap advanced settings',
    'Infomap': 'Infomap advanced settings',
    'htmlTemplate': os.path.join(resourcesDir, 'index.html'),
    'jsTemplate': os.path.join(resourcesDir, 'main.min.js')
}

##############################
# Add variables to namespace #
##############################
for variable,value in global_dict.items():
    setattr(namespace, variable, value)
