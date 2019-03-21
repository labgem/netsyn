#!/usr/bin/env python3

##########
# Import #
##########
import argparse
import json
import logging
import os
import common
import igraph as ig
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

def nodesAttributesInitialization(graph, headersMD=None):
    graph.vs['name'] = []
    graph.vs['id'] = []
    graph.vs['target_idx'] = []
    #graph.vs['context_idx'] = []
    #graph.vs['organism_idx'] = []
    graph.vs['NetSyn_WalkTrap'] = []
    graph.vs['NetSyn_Louvain'] = []
    graph.vs['NetSyn_Infomap'] = []
    graph.vs['UniProt_AC'] = []
    graph.vs['protein_AC'] = []
    graph.vs['organism_name'] = []
    graph.vs['organism_strain'] = []
    graph.vs['organism_taxon_id'] = []
    for rank in common.global_dict['desired_ranks_lneage'].keys():
        graph.vs['lineage_{}'.format(rank)] = []
    return graph

def edgesAttributesInitialization(graph):
    edges_attributes = ['protein_AC', 'products', 'ec_numbers', 'UniProt_AC', 'gene_names', 'locus_tag']
    graph.es['MMseqs_families'] = []
    for attribute in edges_attributes:
        graph.es['{}_source'.format(attribute)] = []
        graph.es['{}_target'.format(attribute)] = []
    return graph, edges_attributes

def fill_node_information(node_reference, prot_in_synteny, prot_family, edges_attributes):
    for pkey, pvalue in prot_in_synteny.items():
        if pkey in edges_attributes:
            node_reference[pkey].setdefault(prot_family, []).append(pvalue)
    return node_reference

def createFullGraph(allData, headersMD=None):
    graph = ig.Graph()
    graph = nodesAttributesInitialization(graph, headersMD)
    graph, edgesAttributes = edgesAttributesInitialization(graph)

    # COM: Generate nodes in a graph
    for idx, node in enumerate(allData['nodes']):
        graph.add_vertex(node['protein_AC'])
        for nkey, nvalue in node.items():
            # COM: addition of organism information to the node
            if nkey == 'organism_idx':
                for okey, ovalue in allData['organisms'][nvalue].items():
                    if okey not in ['id', 'targets_idx', 'lineage']:
                        graph.vs[idx]['organism_{}'.format(okey)] = ovalue
                    elif okey == 'lineage':
                        for level in allData['organisms'][nvalue]['lineage']:
                            rank = level['rank']
                            scientificName = level['scientificName']
                            graph.vs[idx]['lineage_{}'.format(rank)] = scientificName
            # COM: writting every cluster out of a list
            elif nkey == 'clusterings':
                for ckey, cvalue in node[nkey].items():
                    graph.vs[idx]['NetSyn_{}'.format(ckey)] = cvalue
            # COM: recovery of metadata information //// lists of metadata do not need to be initialized :/ very strange ...
            elif nkey == 'metadata':
                for mkey, mvalue in node[nkey].items():
                    graph.vs[idx]['metadata_{}'.format(mkey)] = mvalue
            # COM: recovery of identifiers and accession numbers
            elif nkey in ['id', 'target_idx', 'UniProt_AC', 'protein_AC']:
                graph.vs[idx][nkey] = nvalue

    # COM: Generate edges in a graph
    for edge in allData['edges']:
        graph.add_edge(edge['source'], edge['target'])
        edge_index = graph.get_eid(edge['source'], edge['target'])
        graph.es[edge_index]['weight'] = edge['weight']

        # COM: generate a dictionnary to get all information of the proteins involved in the synteny relation between two 'targets' (source, target)
        attribute_families = {
            'source': {},
            'target': {}
            }
        for attr in edgesAttributes:
            attribute_families['source'].setdefault(attr, {})
            attribute_families['target'].setdefault(attr, {})
        for prot_idx in edge['proteins_idx']:
            prot_in_synt = allData['proteins'][prot_idx]
            prot_family = prot_in_synt['family']
            for prot_target in prot_in_synt['targets_idx']:
                if prot_target == allData['nodes'][edge['source']]['target_idx']:
                    attribute_families['source'] = fill_node_information(attribute_families['source'], prot_in_synt, prot_family, edgesAttributes)
                    break
                elif prot_target == allData['nodes'][edge['target']]['target_idx']:
                    attribute_families['target'] = fill_node_information(attribute_families['target'], prot_in_synt, prot_family, edgesAttributes)
                    break
                else:
                    print('gros probleme')

        # COM: concatenate all the values with ' ~~ ' as separator between values belonging to the same MMseqs family and ' |-| ' as separator between families
        for node_type, node_information in attribute_families.items(): # node_type = source OR target
            for attr_name, attr in node_information.items(): # attr_name = one after the other value of edgesAttibutes
                complete_value = []
                families = list(attr.keys())
                families.sort()
                for afam in families:
                    glued_content = ' ~~ '.join(attr[afam])
                    complete_value.append(glued_content)
                graph.es[edge_index]['{}_{}'.format(attr_name, node_type)] = ' |-| '.join(complete_value)
        graph.es[edge_index]['MMseqs_families'] = ' |-| '.join(str(fam) for fam in families)
    return graph

def run(nodesFile, edgesFile, organismsFile, proteinsFile, metadataFile, resultDir):
    # Constants
    boxName = common.global_dict['boxName']['DataExport']
    # Outputs
    graphmlOut = common.global_dict['files']['DataExport']['graphML']
    htmlOut = common.global_dict['files']['DataExport']['html']
    #settingsOut = common.global_dict['files']['DataExport']['settings']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    print('')
    logger.info('{} running...'.format(boxName))
    # Process
    nodesContent = common.readJSON(nodesFile)
    if metadataFile:
        metadataContent, headersMD = checkAndGetMetadata(metadataFile)
        nodesContent = insertMetadata(nodesContent, metadataContent, headersMD)

    edgesContent = common.readJSON(edgesFile)
    organismsContent = common.readJSON(organismsFile)
    proteinsContent = common.readJSON(proteinsFile)

    if not os.path.isdir(resultDir):
        os.mkdir(resultDir)

    netsynResult = {
        'nodes': nodesContent,
        'edges': edgesContent,
        'proteins': proteinsContent,
        'organisms': organismsContent
    }
    common.write_json(netsynResult, htmlOut)

    if metadataFile:
        full_graph = createFullGraph(netsynResult, headersMD)
    else:
        full_graph = createFullGraph(netsynResult)

    full_graph.write_graphml(graphmlOut)
    logger.info('End of DataEport')
    logger.info('Last advice for life: enjoy your work !!!')

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                     usage='''DataExport options...''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('--nodesFile', type=str,
                        required=True, help='Path of the nodes file obtained from the SyntenyFinder part')
    group1.add_argument('--edgesFile', type=str,
                        required=True, help='Path of the edges file obtained from the SyntenyFinder part')
    group1.add_argument('--organismsFile', type=str,
                        required=True, help='Path of the organims file from the ClusteringIntoFamilies part')
    group1.add_argument('--proteinsFile', type=str,
                        required=True, help='Path of the proteins file from the ClusteringIntoFamilies part')
    group1.add_argument('--metadataFile', type=str,
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
    run(args.nodesFile, args.edgesFile, args.organismsFile, args.proteinsFile, args.metadataFile, '.')
