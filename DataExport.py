#!/usr/bin/env python3

##########
# Import #
##########
import argparse
import json
import logging
import os
import common
#############
# Functions #
#############
def checkAndGetMetadata(metadataFileName):
    '''
    Valide or unvalide the metadata file.
    Initializes to "default value" the undefined metadata.
    '''
    logger = logging.getLogger('{}.{}'.format(checkAndGetMetadata.__module__, checkAndGetMetadata.__name__))
    sep = '\t'
    error = False
    metadataContent = {}
    with open(metadataFileName, 'r') as file:
        nbLines = 0
        for line in file:
            nbLines += 1
            values = [value.rstrip() for value in line.split(sep)]
            if nbLines == 1:
                headersMD = {}
                nbHeaders = len(values)
                for index, value in enumerate(values):
                    if value in headersMD.values():
                        logger.error('{}: Duplicated column.'.format(value))
                        error = True
                    elif value == '':
                        logger.error('Empty field at line of headers.')
                        error = True
                    else:
                        headersMD[index] = value
                        if value == 'accession_type':
                            accessionTypeIndex = index
                for mandatorycolumn in common.global_dict['metadataMadatoryColumn']:
                    if not mandatorycolumn in headersMD.values():
                        logger.error('Missing column, {} column is mandatory.'.format(mandatorycolumn))
                        error = True
                if error:
                    break
            else:
                if not len(values) == nbHeaders:
                    logger.error('At line {}: the number of columns is not egal than the number of headers.'.format(nbLines))
                    error = True
                    continue
                metadata = {}
                for index, value in enumerate(values):
                    if index == accessionTypeIndex and value not in common.global_dict['metadataAccessionAuthorized']:
                        logger.error('At line {}: the "accession_type" must be egal to {}.'.format(nbLines,' or '.join(common.global_dict['metadataAccessionAuthorized'])))
                        error = True
                    else:
                        if headersMD[index] == 'accession_type':
                            continue
                        elif headersMD[index] == 'accession':
                            currentAccession = value
                        else:
                            metadata[headersMD[index]] = common.global_dict['defaultValue'] if value == '' else value
                metadataContent[currentAccession] = metadata
    if error:
        logger.error('Improper metadata file.')
        exit(1)
    return metadataContent, headersMD

def insertMetadata(nodesContent, metadataContent, headersMD):
    for node in nodesContent:
        if node[common.global_dict['inputIheader']] not in metadataContent.keys() and node[common.global_dict['proteinACHeader']] not in metadataContent.keys():
            node['metadata'] = {}
            for header in headersMD.values():
                if header != 'accession_type' and header != 'accession':
                    node['metadata'][header] = common.global_dict['defaultValue']
        else:
            if node[common.global_dict['inputIheader']] in metadataContent.keys():
                accession = node[common.global_dict['inputIheader']]
            elif node[common.global_dict['proteinACHeader']] in metadataContent.keys():
                accession = node[common.global_dict['proteinACHeader']]
            node['metadata'] = metadataContent[accession]
    return nodesContent

def run(nodesFile, edgesFile, organismsFile, proteinsFile, metadataFile, resultDir):
    # Constants
    boxName = common.global_dict['boxName']['DataExport']
    # Outputs
    graphML = common.global_dict['files']['DataExport']['graphML']
    htmlOut = common.global_dict['files']['DataExport']['html']
    #settingsOut = common.global_dict['files']['DataExport']['settings']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    print('')
    logger.info('{} running...'.format(boxName))
    # Process
    metadataContent, headersMD = checkAndGetMetadata(metadataFile) if metadataFile else {}
    nodesContent = common.readJSON(nodesFile)
    nodesContent = insertMetadata(nodesContent, metadataContent, headersMD)

    edgesContent = common.readJSON(edgesFile)
    organismsContent = common.readJSON(organismsFile)
    proteinsContent = common.readJSON(proteinsFile)

    # print('> nodesContent')
    # for i in nodesContent:
    #     print(i)
    # print('------------------')

    # print('> edgesContent')
    # for i in edgesContent:
    #     print(i)
    # print('------------------')

    # print('> organismsContent')
    # for i in organismsContent:
    #     print(i)
    # print('------------------')

    # print('> proteinsContent')
    # for i in proteinsContent:
    #     print(i)
    # print('------------------')

    if not os.path.isdir(resultDir):
        os.mkdir(resultDir)

    netsynResult = {
        'nodes': nodesContent,
        'edges': edgesContent,
        'proteins': proteinsContent,
        'organisms': organismsContent
    }
    common.write_json(netsynResult, htmlOut)

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                     usage='''DataExport options...''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('--Nodes', type=str,
                        required=True, help='Path of the nodes file obtained from the SyntenyFinder part')
    group1.add_argument('--Edges', type=str,
                        required=True, help='Path of the edges file obtained from the SyntenyFinder part')
    group1.add_argument('--Taxonomy', type=str,
                        required=True, help='Path of the taxonomyLineage file from the ClusteringIntoFamilies part')
    group1.add_argument('--ContigsInfo', type=str,
                        required=True, help='Path of the contigs file from the ClusteringIntoFamilies part')
    group1.add_argument('--MetaData', type=str,
                        required=True, help='Path of the metadata file provided by the user')
    group1.add_argument('--OutputName', type=str, required=True,
                        help='Output name files.')

    group2 = parser.add_argument_group('General settings')
    group2.add_argument('--UniProtACList', type=str,
                        help='UniProt accession list input(cf: wiki).')
    group2.add_argument('--CorrespondingFile', type=str,
                        help='Input file of corresponding between: protein_AC/nucleic_AC/nucleic_File_Path (cf: wiki).')

    group3 = parser.add_argument_group('logger')
    group3.add_argument( '--log_level',
                         type = str,
                         nargs = '?',
                         default = 'INFO',
                         help = 'log level',
                         choices = ['ERROR', 'error', 'WARNING', 'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                         required = False )
    group3.add_argument( '--log_file',
                         type = str,
                         nargs = '?',
                         help = 'log file (use the stderr by default)',
                         required = False )
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
    #############
    # Constants #
    #############
    boxName = common.global_dict['boxName']['DataExport']
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('graphML', '{}_Results.graphML'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('html', '{}_Results.html'.format(args.OutputName))
    #######
    # Run #
    #######
    run(args.Nodes, args.Edges, args.Taxonomy, args.ContigsInfo, args.MetaData, '.')
