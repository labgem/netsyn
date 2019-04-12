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
                        logger.error('{}: Duplicated column'.format(value))
                        error = True
                    elif value == '':
                        logger.error('Header column without name')
                        error = True
                    else:
                        headersMD[index] = value
                        if value == 'accession_type':
                            accessionTypeIndex = index
                for mandatorycolumn in common.global_dict['metadataMadatoryColumn']:
                    if not mandatorycolumn in headersMD.values():
                        logger.error('Missing column, {} column is mandatory'.format(mandatorycolumn))
                        error = True
                if error:
                    break
            else:
                if not len(values) == nbHeaders:
                    logger.error('Line {}: a field is empty. Please fulfill the field with a {} if no value'.format(nbLines, common.global_dict['defaultValue']))
                    error = True
                    continue
                metadata = {}
                for index, value in enumerate(values):
                    if index == accessionTypeIndex and value not in common.global_dict['metadataAccessionAuthorized']:
                        logger.error('Line {}: the "accession_type" value must be: {}'.format(nbLines,' or '.join(common.global_dict['metadataAccessionAuthorized'])))
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
        logger.error('Improper metadata file')
        exit(1)
    return metadataContent, headersMD

def insertMetadata(nodesContent, metadataContent, headersMD):
    '''
    Insertes the metada into nodes content.
    '''
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

def getTaxID(rank, organismIndex, organismsContent):
    '''
    Get the taxonomic id for the rank of the organism.
    '''
    logger = logging.getLogger('{}.{}'.format(getTaxID.__module__, getTaxID.__name__))
    for level in organismsContent[organismIndex]['lineage']:
        if level['rank'] == rank:
            return level['tax_id']
    logger.critical('Tax_id corresponding to the {} rank not found for the {} {} organism (index: {})'.format(rank,
                                                                                                               organismsContent[organismIndex]['name'],
                                                                                                               organismsContent[organismIndex]['strain'],
                                                                                                               organismIndex))
    exit(1)

def getNodesToMerge(nodes, criterionType, criterion, clusteringMethod, organismsContent=None):
    '''
    Determines which nodes to merge.
    '''
    nodesPerCriterion = {}
    for index, node in enumerate(nodes):
        if criterionType == 'taxonomic':
            node.setdefault(criterionType, {}).setdefault(criterion, getTaxID(criterion, node['organism_idx'], organismsContent))
        if not node[criterionType][criterion] == common.global_dict['defaultValue']:
            node['node_idx'] = index
            nodesPerCriterion.setdefault(node['clusterings'][clusteringMethod], {}).setdefault(node[criterionType][criterion], []).append(node)
        if criterionType == 'taxonomic':
            del node[criterionType]
    nodesToMerge = []
    for cluster in nodesPerCriterion:
        for nodes in nodesPerCriterion[cluster]:
            if len(nodesPerCriterion[cluster][nodes]) > 1:
                nodesToMerge.append(nodesPerCriterion[cluster][nodes])
    return nodesToMerge

def nodes_organismsMerging(nodesToMerge, organismsContent, clusteringMethod):
    '''
    Merges nodes content and updates organisms content.
    '''
    logger = logging.getLogger('{}.{}'.format(nodes_organismsMerging.__module__, nodes_organismsMerging.__name__))
    newNodes = []
    organismIdMax = max([organism['id'] for organism in organismsContent])
    for nodes in nodesToMerge:
        newNode = {}
        organismToMerge = []
        oldNodes = []
        metadata = {}
        clusterings = {}
        clusterings[clusteringMethod] = None
        for node in nodes:
            oldNodes.append(node['node_idx'])
            newNode.setdefault('id', []).append(str(node['id']))
            newNode.setdefault('target_idx', []).append(str(node['target_idx']))
            newNode.setdefault(common.global_dict['inputIheader'], []).append(node[common.global_dict['inputIheader']])
            newNode.setdefault(common.global_dict['proteinACHeader'], []).append(node[common.global_dict['proteinACHeader']])
            if not clusterings[clusteringMethod]:
                clusterings[clusteringMethod] = node['clusterings'][clusteringMethod]
            elif clusterings[clusteringMethod] != node['clusterings'][clusteringMethod]:
                logger.critical('Nodes to merge must belong to the same cluster')
            if 'metadata' in node.keys():
                for key, value in node['metadata'].items():
                    if key not in metadata:
                        metadata[key] = []
                    if value not in metadata[key]:
                        metadata[key].append(value)
            if organismsContent[node['organism_idx']] not in organismToMerge:
                organismToMerge.append(organismsContent[node['organism_idx']])
        for key, toMerge in newNode.items():
            newNode[key] = ', '.join(toMerge)
        newNode['oldNodes'] = oldNodes
        if metadata != {}:
            for key, toMerge in metadata.items():
                metadata[key] = ', '.join(sorted(toMerge))
            newNode['metadata'] = metadata
        newNode['clusterings'] = clusterings
        newNode['Size'] = len(nodes)

        metaLineage = {}
        organismIdMax += 1
        newOrganism = {'id': organismIdMax}
        for organism in organismToMerge:
            newOrganism.setdefault('strain', []).append(organism['strain'])
            newOrganism.setdefault('name', []).append(organism['name'])
            newOrganism.setdefault('taxon_id', []).append(organism['taxon_id'])
            newOrganism.setdefault('lineage', [])
            for level in organism['lineage']:
                if level['level'] not in metaLineage:
                    metaLineage[level['level']] = {
                        'rank': level['rank'],
                        'scientificName': [],
                        'tax_id': []
                    }
                if level['scientificName'] not in metaLineage[level['level']]['scientificName']:
                    metaLineage[level['level']]['scientificName'].append(level['scientificName'])
                if str(level['tax_id']) not in metaLineage[level['level']]['tax_id']:
                    metaLineage[level['level']]['tax_id'].append(str(level['tax_id']))
        for level, values in metaLineage.items():
            newOrganism['lineage'].append({
                'rank': values['rank'],
                'scientificName': ', '.join(values['scientificName']),
                'tax_id': ', '.join(values['tax_id']),
                'level': level
            })
        for key in ['strain', 'name', 'taxon_id']:
            newOrganism[key] = ', '.join(newOrganism[key])
        organismsContent.append(newOrganism)
        newNode['organism_id'] = organismIdMax
        newNode['organism_idx'] = (len(organismsContent)-1)
        newNodes.append(newNode)
    return newNodes, organismsContent

def nodes_edgesUpdating(nodesToAdd, nodesContent, edgesContent, clusteringMethod):
    '''
    Addes the nodes into nodes list and updates edges content.
    '''
    nodesIndexes = [i for i in range(len(nodesContent))]
    nodesToDelIndexes = []
    for nodeToAdd in nodesToAdd:
        for nodeToDel in nodeToAdd['oldNodes']:
            nodesToDelIndexes.append(nodeToDel)
    decrement = 0
    for i in nodesIndexes:
        if i in nodesToDelIndexes:
            nodesIndexes[i] = -1
            decrement += 1
        else:
            nodesIndexes[i] -= decrement
    for nodeToDel in sorted(nodesToDelIndexes, reverse=True):
        del nodesContent[nodeToDel]
    for nodeToAdd in nodesToAdd:
        nodesToDel = nodeToAdd['oldNodes']
        del nodeToAdd['oldNodes']
        newNodeIndex = len(nodesContent)
        nodesContent.append(nodeToAdd)
        for nodeToDel in nodesToDel:
            nodesIndexes[nodeToDel] = newNodeIndex
    for node in nodesContent:
        clusterID = node['clusterings'][clusteringMethod]
        node['clusterings'] = {}
        node['clusterings'][clusteringMethod] = clusterID

    edgesToDel = []
    for index,edge in enumerate(edgesContent):
        edge['source'] = nodesIndexes[int(edge['source'])]
        edge['target'] = nodesIndexes[int(edge['target'])]
        if edge['source'] == edge['target']:
            edgesToDel.append(index)
    for edgeIndex in sorted(edgesToDel, reverse=True):
        del edgesContent[edgeIndex]
    return nodesContent, edgesContent

def fill_node_information(node_reference, prot_in_synteny, prot_family, edges_attributes):
    for pkey, pvalue in prot_in_synteny.items():
        if pkey in edges_attributes:
            node_reference[pkey].setdefault(prot_family, []).append(pvalue)
    return node_reference

def get_neighbour(proteins_idx, edgesAttributes, attribute_families, allData):
    for prot_idx in proteins_idx:
        prot_in_synt = allData['proteins'][int(prot_idx)]
        prot_family = prot_in_synt['family']
        attribute_families = fill_node_information(attribute_families, prot_in_synt, prot_family, edgesAttributes)
    return attribute_families

def createFullGraph(allData, clusteringMethod, headersMD=None):
    graph = ig.Graph()
    edgesAttributes = common.global_dict['graphML_edgesAttributes']

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
                if clusteringMethod:
                    graph.vs[idx]['NetSyn_{}'.format(clusteringMethod)] = node[nkey][clusteringMethod]
                else:
                    for ckey, cvalue in node[nkey].items():
                        graph.vs[idx]['NetSyn_{}'.format(ckey)] = cvalue
            # COM: recovery of metadata information //// lists of metadata do not need to be initialized :/ very strange ...
            elif nkey == 'metadata':
                for mkey, mvalue in node[nkey].items():
                    graph.vs[idx]['metadata_{}'.format(mkey)] = mvalue
            # COM: recovery of identifiers and accession numbers
            elif nkey in ['id', 'target_idx', 'UniProt_AC', 'protein_AC','Size']:
                graph.vs[idx][nkey] = nvalue

    # COM: Generate edges in a graph
    for edge in allData['edges']:
        graph.add_edge(int(edge['source']), int(edge['target']))
        edge_index = graph.get_eid(int(edge['source']), int(edge['target']))
        graph.es[edge_index]['weight'] = edge['weight']

        # COM: generate a dictionnary to get all information of the proteins involved in the synteny relation between two 'targets' (source, target)
        attribute_families = {
            'source': {},
            'target': {}
            }
        for attr in edgesAttributes:
            attribute_families['source'].setdefault(attr, {})
            attribute_families['target'].setdefault(attr, {})

        attribute_families['source'] = get_neighbour(edge['proteins_idx_source'], edgesAttributes, attribute_families['source'], allData)
        attribute_families['target'] = get_neighbour(edge['proteins_idx_target'], edgesAttributes, attribute_families['target'], allData)

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

def run(nodesFile, edgesFile, organismsFile, proteinsFile, metadataFile, redundancyRemovalLabel, redundancyRemovalTaxonomy, clusteringMethod):
    # Constants
    boxName = common.global_dict['boxName']['DataExport']
    # Outputs
    graphmlOut = common.global_dict['files']['DataExport']['graphML']
    htmlOut = common.global_dict['files']['DataExport']['html']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    reportingMessages = []
    print('')
    logger.info('{} running...'.format(boxName))
    # Process
    nodesContent = common.readJSON(nodesFile)
    if metadataFile:
        if not common.checkFilledFile(metadataFile):
            metadataContent, headersMD = checkAndGetMetadata(metadataFile)
            nodesContent = insertMetadata(nodesContent, metadataContent, headersMD)
        else:
            logger.error('Please make sure that {} file is in the appropriate repertory'.format(metadataFile))
            exit(1)

    edgesContent = common.readJSON(edgesFile)
    organismsContent = common.readJSON(organismsFile)
    proteinsContent = common.readJSON(proteinsFile)

    if redundancyRemovalLabel or redundancyRemovalTaxonomy:
        logger.info('Redundancy removal processing...')
        nodeOldsNumber = len(nodesContent)
        if redundancyRemovalLabel:
            if redundancyRemovalLabel not in headersMD.values():
                logger.error('The label {} is not in the metadata column headers: {}'.format(redundancyRemovalLabel, ', '.join(headersMD.values())))
                exit(1)
            reportingMessages.append('Redundancy removal settings used: label "{}", clustering method "{}"'.format(
                redundancyRemovalLabel, clusteringMethod
            ))
            nodesToMerge = getNodesToMerge(nodesContent, 'metadata', redundancyRemovalLabel, clusteringMethod)
        elif redundancyRemovalTaxonomy:
            reportingMessages.append('Redundancy removal settings used: taxonomic rank "{}", clustering method "{}"'.format(
                redundancyRemovalTaxonomy, clusteringMethod
            ))
            nodesToMerge = getNodesToMerge(nodesContent, 'taxonomic', redundancyRemovalTaxonomy, clusteringMethod, organismsContent)
        newNodes, organismsContent = nodes_organismsMerging(nodesToMerge, organismsContent, clusteringMethod)
        nodesContent, edgesContent = nodes_edgesUpdating(newNodes, nodesContent, edgesContent, clusteringMethod)
        reportingMessages.append('Nodes merged after redundancy removal step: {}/{}'.format(len(nodesContent), nodeOldsNumber))

    netsynResult = {
        'nodes': nodesContent,
        'edges': edgesContent,
        'proteins': proteinsContent,
        'organisms': organismsContent
    }

    #####
    # for edge in edgesContent:
    #     logger.debug('{} > {}'.format(edge["source"],edge["target"]))
    # logger.debug('')
    # print('\n\n\n')
    # for n in nodesContent:
    #     print('{}: '.format(n))
    #     # print('{}: '.format(n['target_idx']))
    #     # for l in organismsContent[n['organism_idx']]['lineage']:
    #         # print('    {}'.format(l))
    # print()
    # for e in edgesContent:
    #     s = nodesContent[e['source']]['target_idx']
    #     t = nodesContent[e['target']]['target_idx']
    #     print('{}-{}:{}'.format(s,t,e))
    #####

    logger.info('Output graph building...')
    if metadataFile:
        full_graph = createFullGraph(netsynResult, clusteringMethod, headersMD)
    else:
        full_graph = createFullGraph(netsynResult, clusteringMethod)

    dataDirectoryProcess = os.path.join(common.global_dict['dataDirectory'], boxName)
    if not os.path.isdir(dataDirectoryProcess):
        os.mkdir(dataDirectoryProcess)

    common.write_json(netsynResult, htmlOut)
    full_graph.write_graphml(graphmlOut)
    logger.info('{} completed!'.format(boxName))
    reportingMessages.append('Aller faut trouver un truc a mettre...'.format())
    common.reportingFormat(logger, boxName, reportingMessages)

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                     usage='''DataExport.py -nf <nodesFile> -ef <edgesFile> -of <organismsFile> -pf <proteinsFile> -o <OutputName>\n\
\t\t[-mf <metadataFile>]\n\t\t[-rrl <RedundancyRemovalLabel> | -rrt <RedundancyRemovalTaxonomy> -cm <ClusteringMethod>]''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-nf', '--nodesFile', type=str,
                        required=True, help='Path of the nodes file obtained from the SyntenyFinder part')
    group1.add_argument('-ef', '--edgesFile', type=str,
                        required=True, help='Path of the edges file obtained from the SyntenyFinder part')
    group1.add_argument('-of', '--organismsFile', type=str,
                        required=True, help='Path of the organims file from the ClusteringIntoFamilies part')
    group1.add_argument('-pf', '--proteinsFile', type=str,
                        required=True, help='Path of the proteins file from the ClusteringIntoFamilies part')
    group1.add_argument('-o', '--OutputName', type=str, required=True,
                        help='Output name files')
    group1.add_argument('-mf', '--metadataFile', type=str,
                        required=False, help='Path of the metadata file provided by the user')

    group2 = parser.add_argument_group('Redundancy Removal settings')
    group2.add_argument('-cm', '--ClusteringMethod', type=str,
                        choices=['MCL','Infomap','Louvain','WalkTrap'],
                        default=None,
                        help='Clustering method choose in : MCL (small graph), Infomap (medium graph), Louvain (medium graph) or WalkTrap (big  graph).\nDefault value: MCL')
    group2.add_argument('-rrl', '--RedundancyRemovalLabel', type=str,
                        help='Label of the metadata column on which the redundancy will be computed (Incompatible with --RedundancyRemovalTaxonomy option.)')
    group2.add_argument('-rrt', '--RedundancyRemovalTaxonomy', type=str, choices=common.global_dict['desired_ranks_lineage'].keys(),
                        help='Taxonomic rank on which the redundancy will be computed. (Incompatible with --RedundancyRemovalLabel option)')

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
    return parser.parse_args(), parser

if __name__ == '__main__':
    import argparse
    ######################
    # Parse command line #
    ######################
    args, parser = argumentsParser()
    if args.RedundancyRemovalLabel and args.RedundancyRemovalTaxonomy:
        parser.error('RedundancyRemovalLabel and RedundancyRemovalTaxonomy are incompatible options. Please choose one of two options')
    if (args.RedundancyRemovalLabel or args.RedundancyRemovalTaxonomy) and not args.ClusteringMethod:
        parser.error('A clustering method must be provided since a redundancy removal option (-rrl, -rrt) is selected')
    if args.ClusteringMethod and not (args.RedundancyRemovalLabel or args.RedundancyRemovalTaxonomy):
        parser.error('The selection of a clustering algorithm is available only when redundancy removal is enabled (see DataExport.py -h/--help)')
    if args.RedundancyRemovalLabel and not args.metadataFile:
        parser.error('Please specify the --metadataFile option')
    ##########
    # Logger #
    ##########
    common.parametersLogger(args)
    #############
    # Constants #
    #############
    common.global_dict['dataDirectory'] = '.'
    boxName = common.global_dict['boxName']['DataExport']
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('graphML', '{}_Results.graphML'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('html', '{}_Results.html'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('report', '{}_{}_report.txt'.format(args.OutputName, boxName))
    #######
    # Run #
    #######
    run(args.nodesFile,
        args.edgesFile,
        args.organismsFile,
        args.proteinsFile,
        args.metadataFile,
        args.RedundancyRemovalLabel,
        args.RedundancyRemovalTaxonomy,
        args.ClusteringMethod)
