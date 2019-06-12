#!/usr/bin/env python3

##########
# Import #
##########
import argparse
import logging
import os
import igraph as ig
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
            newNode.setdefault('protein_idx', []).append(str(node['protein_idx']))
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
    for index, edge in enumerate(edgesContent):
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
            elif nkey in ['id', 'target_idx', 'UniProt_AC', 'protein_AC', 'Size']:
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
            # COM: iGraph reorders the source and edge nodes, so the assigment of attributes has to be also inverted
            if int(edge['source']) > int(edge['target']):
                if node_type == 'source':
                    node_type = 'target'
                else:
                    node_type = 'source'
            for attr_name, attr in node_information.items(): # attr_name = one after the other value of edgesAttibutes
                complete_value = []
                families = list(attr.keys())
                if common.global_dict['defaultValue'] in families:
                    families.remove(common.global_dict['defaultValue'])
                    families.sort()
                    families.append(common.global_dict['defaultValue'])
                else:
                    families.sort()
                for afam in families:
                    glued_content = ' ~~ '.join(attr[afam])
                    complete_value.append(glued_content)
                graph.es[edge_index]['{}_{}'.format(attr_name, node_type)] = ' |-| '.join(complete_value)
        graph.es[edge_index]['MMseqs_families'] = ' |-| '.join(str(fam) for fam in families)
    return graph

def get_families_of_targets(nodesContent, proteinsContent):
    families_of_targetsIdx = {}
    for node_idx, node in enumerate(nodesContent):
        for index in str(node['protein_idx']).split(', '):
            family = proteinsContent[int(index)]['family']
            families_of_targetsIdx.setdefault(family, []).append(node_idx)
    return families_of_targetsIdx

def get_clusterID_family_content(node_dict, proteins_indexes, proteinsContent):
    families = []
    for prot_idx in proteins_indexes:
        family = proteinsContent[int(prot_idx)]['family']
        families.append(family)
        node_dict.setdefault(family, {}).setdefault('protein_AC', []).append(proteinsContent[int(prot_idx)]['protein_AC'])
        node_dict[family].setdefault('locus_tags', []).append(proteinsContent[int(prot_idx)]['locus_tag'])
        node_dict[family].setdefault('ec_numbers', []).append(proteinsContent[int(prot_idx)]['ec_numbers'])
        node_dict[family].setdefault('gene_names', []).append(proteinsContent[int(prot_idx)]['gene_names'])
        node_dict[family].setdefault('products', []).append(proteinsContent[int(prot_idx)]['products'])
    node_dict.setdefault('families', families)
    return node_dict, families

def get_strain_and_species(node_idx, taxo_level, organismsContent, nodesContent):
    node_organism = {}
    # COM: strain treatment
    strain_value = organismsContent[nodesContent[node_idx]['organism_idx']]['strain']
    if strain_value == 'NA':
        strain_value = ''
    # COM: species treatment
    for taxonomic_level in organismsContent[nodesContent[node_idx]['organism_idx']]['lineage']:
        if taxonomic_level['rank'] == taxo_level:
            species_value = taxonomic_level['scientificName']
    if species_value == 'NA':
        if organismsContent[nodesContent[node_idx]['organism_idx']]['name'] != 'NA':
            species_value = organismsContent[nodesContent[node_idx]['organism_idx']]['name']
        elif strain_value != '':
            species_value = 'unknown'
        else:
            species_value = 'NA'
    # COM: organism treatment
    organism_value = organismsContent[nodesContent[node_idx]['organism_idx']]['name']
    if strain_value not in organism_value:
        organism_value = ' '.join([organism_value, strain_value])
    # COM: data insertion in node_organism
    node_organism.setdefault('strains', []).append(strain_value)
    node_organism.setdefault('species', []).append(species_value)
    node_organism.setdefault('organisms', []).append(organism_value)
    return node_organism

def get_metadata_values(dico, node_idx, nodesContent):
    for key_md, value_md in nodesContent[node_idx]['metadata'].items():
        dico.setdefault('metadata', {}).setdefault(key_md, []).append(value_md)
    return dico

def get_non_redundant_prots_info(storage_family_dict, cluster_family_dict):
    for prot_idx, protAC in enumerate(storage_family_dict['protein_AC']):
        if protAC not in cluster_family_dict.setdefault('protein_AC', []):
            for key, value in storage_family_dict.items():
                cluster_family_dict.setdefault(key, []).append(value[prot_idx])
    return cluster_family_dict

def add_metadataHeaders(headers_to_print, headersMD):
    metadataHeaders = sorted([header
                              for header in headersMD.values()
                              if header not in common.global_dict['metadataMadatoryColumn']]
                             )
    return metadataHeaders

def get_metadata_info(familyContent, md_dict):
    for key, value in md_dict.items():
        md_val = []
        for elt in value:
            md_val.extend(elt.split(', '))
        familyContent.setdefault('metadata', {}).setdefault(key, []).extend(md_val)
    return familyContent

def get_sorted_nr_data(specific_list, counter=False):
    if counter:
        lower_case_list = [name.lower() for name in specific_list]
        convert_to_dict = {name: lower_case_list.count(name) for name in set(lower_case_list)}
        if 'na' in convert_to_dict.keys() and len(convert_to_dict) == 1:
            sorted_nr_list = ['NA']
        else:
            if 'na' in convert_to_dict.keys():
                del convert_to_dict['na']
            sorted_keys_list = sorted(convert_to_dict.keys(), key=lambda x: x.lower())
            sorted_nr_list = ['{} [{}]'.format(name, convert_to_dict[name]) for name in sorted_keys_list]
    else:
        sorted_nr_list = sorted(set([name for name in specific_list if name != 'NA']), key=lambda x: x.lower())#, key=lambda x: familyContent['gene_names'].index(x))
    string = ' / '.join(sorted_nr_list)
    return string

def add_metadata_info(line, metadataHeaders, headers_to_print, familyContent):
    for md_head in metadataHeaders:
        line[headers_to_print.index(md_head)] = ' / '.join(set(familyContent['metadata'][md_head]))
    return line

def skip_duplicates(iterable, key=lambda x: x):
    ''' remove duplicates from a list keeping the order of the elements
    Use a generator
    '''
    # on va mettre l’empreinte unique de chaque élément dans ce set
    fingerprints = set()
    for x in iterable:
        # chaque élement voit son emprunte calculée. Par défaut l’empreinte
        # est l'élément lui même, ce qui fait qu'il n'y a pas besoin de
        # spécifier 'key' pour des primitives comme les ints ou les strings.
        fingerprint = key(x)
        # On vérifie que l'empreinte est dans la liste des empreintes  des
        # éléments précédents. Si ce n'est pas le cas, on yield l'élément, et on
        # rajoute sont empreinte ans la liste de ceux trouvés, donc il ne sera
        # pas yieldé si on ne le yieldera pas une seconde fois si on le
        # rencontre à nouveau
        if fingerprint not in fingerprints:
            yield x
            fingerprints.add(fingerprint)

def run(nodesFile, edgesFile, organismsFile, proteinsFile, metadataFile, redundancyRemovalLabel, redundancyRemovalTaxonomy, clusteringMethod):
    # Constants
    boxName = common.global_dict['boxName']['DataExport']
    # Outputs
    graphmlOut = common.global_dict['files']['DataExport']['graphML']
    htmlOut = common.global_dict['files']['DataExport']['html']
    synthesisDirectory = common.global_dict['synthesisDataExport']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    reportingMessages = []
    print('')
    logger.info('{} running...'.format(boxName))
    # Process
    nodesContent = common.readJSON(nodesFile, common.getNodesListStepschema())
    edgesContent = common.readJSON(edgesFile, common.getEdgesListStepschema())
    organismsContent = common.readJSON(organismsFile, common.getOrganismsTaxonomyStepschema())
    proteinsContent = common.readJSON(proteinsFile, common.getProteinsSyntenyStepSchema())
    if metadataFile:
        if not common.checkFilledFile(metadataFile):
            metadataContent, headersMD = checkAndGetMetadata(metadataFile)
            nodesContent = insertMetadata(nodesContent, metadataContent, headersMD)
        else:
            logger.error('Please make sure that {} file is in the appropriate repertory'.format(metadataFile))
            exit(1)

    if redundancyRemovalLabel or redundancyRemovalTaxonomy:
        logger.info('Redundancy removal processing...')
        nodeOldsNumber = len(nodesContent)
        if redundancyRemovalLabel:
            if redundancyRemovalLabel not in headersMD.values():
                logger.error('The label {} is not in the metadata column headers: {}'.format(redundancyRemovalLabel, ', '.join(set(headersMD.values())-set(common.global_dict['metadataMadatoryColumn']))))
                exit(1)
            elif redundancyRemovalLabel in common.global_dict['metadataMadatoryColumn']:
                logger.error('The label {} is not considered as a metadata label. Please choose another label from: {}'.format(redundancyRemovalLabel, ', '.join(set(headersMD.values())-set(common.global_dict['metadataMadatoryColumn']))))
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

    logger.info('Output graph building...')
    if metadataFile:
        full_graph = createFullGraph(netsynResult, clusteringMethod, headersMD)
    else:
        full_graph = createFullGraph(netsynResult, clusteringMethod)

    dataDirectoryProcess = os.path.join(common.global_dict['dataDirectory'], boxName)
    if not os.path.isdir(dataDirectoryProcess):
        os.mkdir(dataDirectoryProcess)
    if not os.path.isdir(synthesisDirectory):
        os.mkdir(synthesisDirectory)

    common.write_json(netsynResult, htmlOut)
    full_graph.write_graphml(graphmlOut)

    ### COM: Storage of information to write ClusteringSynthesis files
    families_of_targetsIdx = get_families_of_targets(nodesContent, proteinsContent)
    intra_cluster = {}
    inter_cluster = {}
    for edge in edgesContent:
        # COM: get the clusterID for each clustering method of source and target nodes
        source_clusterings = nodesContent[int(edge['source'])]['clusterings']
        target_clusterings = nodesContent[int(edge['target'])]['clusterings']
        for cm, clusterID in source_clusterings.items():
            source_dict = {}
            target_dict = {}
            source_dict, source_families = get_clusterID_family_content(source_dict, edge['proteins_idx_source'], proteinsContent)
            target_dict, target_families = get_clusterID_family_content(target_dict, edge['proteins_idx_target'], proteinsContent)
            source_orgs = get_strain_and_species(int(edge['source']), 'species', organismsContent, nodesContent)
            target_orgs = get_strain_and_species(int(edge['target']), 'species', organismsContent, nodesContent)
            if metadataFile:
                source_dict = get_metadata_values(source_dict, int(edge['source']), nodesContent)
                target_dict = get_metadata_values(target_dict, int(edge['target']), nodesContent)
            if clusterID == target_clusterings[cm]:
                if set(source_families)-set(target_families) == set():
                    for family in set(source_families):
                        intra_cluster.setdefault(cm, {}).setdefault(clusterID, {}).setdefault(family, {}).setdefault('nodes_indexes', [])

                        if (int(edge['source']) not in intra_cluster[cm][clusterID][family]['nodes_indexes']) and (int(edge['target']) not in intra_cluster[cm][clusterID][family]['nodes_indexes']):
                            intra_cluster[cm][clusterID][family]['nodes_indexes'].extend([int(edge['source']), int(edge['target'])])
                            intra_cluster[cm][clusterID][family] = get_non_redundant_prots_info(source_dict[family], intra_cluster[cm][clusterID][family])
                            intra_cluster[cm][clusterID][family] = get_non_redundant_prots_info(target_dict[family], intra_cluster[cm][clusterID][family])
                            intra_cluster[cm][clusterID][family].setdefault('strains', []).extend(source_orgs['strains'])
                            intra_cluster[cm][clusterID][family].setdefault('species', []).extend(source_orgs['species'])
                            intra_cluster[cm][clusterID][family].setdefault('organisms', []).extend(source_orgs['organisms'])
                            intra_cluster[cm][clusterID][family].setdefault('strains', []).extend(target_orgs['strains'])
                            intra_cluster[cm][clusterID][family].setdefault('species', []).extend(target_orgs['species'])
                            intra_cluster[cm][clusterID][family].setdefault('organisms', []).extend(target_orgs['organisms'])

                        elif int(edge['source']) not in intra_cluster[cm][clusterID][family]['nodes_indexes']:
                            intra_cluster[cm][clusterID][family]['nodes_indexes'].append(int(edge['source']))
                            intra_cluster[cm][clusterID][family] = get_non_redundant_prots_info(source_dict[family], intra_cluster[cm][clusterID][family])
                            intra_cluster[cm][clusterID][family].setdefault('strains', []).extend(source_orgs['strains'])
                            intra_cluster[cm][clusterID][family].setdefault('species', []).extend(source_orgs['species'])
                            intra_cluster[cm][clusterID][family].setdefault('organisms', []).extend(source_orgs['organisms'])

                        elif int(edge['target']) not in intra_cluster[cm][clusterID][family]['nodes_indexes']:
                            intra_cluster[cm][clusterID][family]['nodes_indexes'].append(int(edge['target']))
                            intra_cluster[cm][clusterID][family] = get_non_redundant_prots_info(target_dict[family], intra_cluster[cm][clusterID][family])
                            intra_cluster[cm][clusterID][family].setdefault('strains', []).extend(target_orgs['strains'])
                            intra_cluster[cm][clusterID][family].setdefault('species', []).extend(target_orgs['species'])
                            intra_cluster[cm][clusterID][family].setdefault('organisms', []).extend(target_orgs['organisms'])

                        intra_cluster[cm][clusterID][family]['synteny_size'] = intra_cluster[cm][clusterID][family].setdefault('synteny_size', 0) + 1
                        if metadataFile:
                            intra_cluster[cm][clusterID][family] = get_metadata_info(intra_cluster[cm][clusterID][family], source_dict['metadata'])
                            intra_cluster[cm][clusterID][family] = get_metadata_info(intra_cluster[cm][clusterID][family], target_dict['metadata'])
                else:
                    logger.critical('Families ID obtained by proteins for both targets in synteny must be identical. This reveals a development error. Please contact us.')
            else:
                families = source_families
                families.extend(target_families)
                for family in set(families):
                    inter_cluster.setdefault(cm, {}).setdefault(family, {}).setdefault('clusters', []).extend([clusterID, target_clusterings[cm]])
                    inter_cluster[cm][family].setdefault('nodes_indexes', [])
                    if (int(edge['source']) and int(edge['target'])) not in inter_cluster[cm][family]['nodes_indexes']:
                        inter_cluster[cm][family]['nodes_indexes'].extend([int(edge['source']), int(edge['target'])])
                        inter_cluster[cm][family] = get_non_redundant_prots_info(source_dict[family], inter_cluster[cm][family])
                        inter_cluster[cm][family] = get_non_redundant_prots_info(target_dict[family], inter_cluster[cm][family])
                        inter_cluster[cm][family].setdefault('strains', []).extend(source_orgs['strains'])
                        inter_cluster[cm][family].setdefault('species', []).extend(source_orgs['species'])
                        inter_cluster[cm][family].setdefault('organisms', []).extend(source_orgs['organisms'])
                        inter_cluster[cm][family].setdefault('strains', []).extend(target_orgs['strains'])
                        inter_cluster[cm][family].setdefault('species', []).extend(target_orgs['species'])
                        inter_cluster[cm][family].setdefault('organisms', []).extend(target_orgs['organisms'])

                    elif int(edge['source']) not in inter_cluster[cm][family]['nodes_indexes']:
                        inter_cluster[cm][family]['nodes_indexes'].append(int(edge['source']))
                        inter_cluster[cm][family] = get_non_redundant_prots_info(source_dict[family], inter_cluster[cm][family])
                        inter_cluster[cm][family].setdefault('strains', []).extend(source_orgs['strains'])
                        inter_cluster[cm][family].setdefault('species', []).extend(source_orgs['species'])
                        inter_cluster[cm][family].setdefault('organisms', []).extend(source_orgs['organisms'])

                    elif int(edge['target']) not in inter_cluster[cm][family]['nodes_indexes']:
                        inter_cluster[cm][family]['nodes_indexes'].append(int(edge['target']))
                        inter_cluster[cm][family] = get_non_redundant_prots_info(target_dict[family], inter_cluster[cm][family])
                        inter_cluster[cm][family].setdefault('strains', []).extend(target_orgs['strains'])
                        inter_cluster[cm][family].setdefault('species', []).extend(target_orgs['species'])
                        inter_cluster[cm][family].setdefault('organisms', []).extend(target_orgs['organisms'])

                    inter_cluster[cm][family]['synteny_size'] = inter_cluster[cm][family].setdefault('synteny_size', 0) + 1
                    if metadataFile:
                        inter_cluster[cm][family] = get_metadata_info(inter_cluster[cm][family], source_dict['metadata'])
                        inter_cluster[cm][family] = get_metadata_info(inter_cluster[cm][family], target_dict['metadata'])

    # COM: intra cluster process
    headers_to_print = ['NetSyn_ClusterID', 'NetSyn_FamilyID', 'Target_Family', 'Syntenies_Nb', 'Species_Nb', 'Strains_Nb', 'Organisms_Nb', 'Organisms_List', 'Proteins_in_Synteny_Nb', 'Products', 'Gene_Names', 'EC_Numbers', 'Locus_Tags', 'Proteins_AC']
    if metadataFile:
        metadataHeaders = add_metadataHeaders(headers_to_print, headersMD)
        headers_to_print.extend(metadataHeaders)
    for cm, cmContent in intra_cluster.items():
        lines_to_print = []
        # COM: generation of lines containing intra cluster information
        for clusterID, clusterContent in cmContent.items():
            for familyID, familyContent in clusterContent.items():
                line = [''] * len(headers_to_print)
                containTarget = 'Y' if familyID in families_of_targetsIdx and set(familyContent['nodes_indexes']) & set(families_of_targetsIdx[familyID]) else 'N'
                strains_nb = len(set(familyContent['strains']))
                species_nb = len(set(familyContent['species']))
                orgs_list = ' / '.join(sorted(set(familyContent['organisms']), key=lambda x: x.lower()))
                orgs_nb = len(set(familyContent['organisms']))

                line[headers_to_print.index('NetSyn_ClusterID')] = clusterID
                line[headers_to_print.index('NetSyn_FamilyID')] = familyID
                line[headers_to_print.index('Target_Family')] = containTarget
                line[headers_to_print.index('Syntenies_Nb')] = familyContent['synteny_size']
                line[headers_to_print.index('Species_Nb')] = species_nb
                line[headers_to_print.index('Strains_Nb')] = strains_nb
                line[headers_to_print.index('Organisms_Nb')] = orgs_nb
                line[headers_to_print.index('Organisms_List')] = orgs_list
                line[headers_to_print.index('Proteins_in_Synteny_Nb')] = len(set(familyContent['protein_AC'])) # use of set() no longer useful
                line[headers_to_print.index('Products')] = get_sorted_nr_data(familyContent['products'], counter=True)
                line[headers_to_print.index('Gene_Names')] = get_sorted_nr_data(familyContent['gene_names'], counter=True)
                line[headers_to_print.index('EC_Numbers')] = get_sorted_nr_data(familyContent['ec_numbers'], counter=True)
                line[headers_to_print.index('Locus_Tags')] = get_sorted_nr_data(familyContent['locus_tags'], counter=False)
                line[headers_to_print.index('Proteins_AC')] = get_sorted_nr_data(familyContent['protein_AC'], counter=False)

                if metadataFile:
                    line = add_metadata_info(line, metadataHeaders, headers_to_print, familyContent)
                lines_to_print.append(line)
        # COM: ordering lines
        order_targetFam = ['Y', 'N']
        order_NSclusters = sorted([(line[headers_to_print.index('NetSyn_ClusterID')],
                                    line[headers_to_print.index('Syntenies_Nb')],
                                    line[headers_to_print.index('Species_Nb')],
                                    line[headers_to_print.index('Strains_Nb')],
                                    line[headers_to_print.index('Proteins_in_Synteny_Nb')])
                                   for line in lines_to_print
                                   if line[headers_to_print.index('Target_Family')] == 'Y'],
                                  key=lambda figures: (-figures[1], -figures[2], -figures[3], -figures[4]))
        order_NSclusters = [figures[0] for figures in order_NSclusters]
        order_NSclusters = list(skip_duplicates(order_NSclusters))

        sorted_lines = sorted(lines_to_print, key=lambda line: (order_NSclusters.index(line[headers_to_print.index('NetSyn_ClusterID')]),
                                                                order_targetFam.index(line[headers_to_print.index('Target_Family')]),
                                                                -line[headers_to_print.index('Syntenies_Nb')],
                                                                -line[headers_to_print.index('Species_Nb')],
                                                                -line[headers_to_print.index('Strains_Nb')],
                                                                -line[headers_to_print.index('Proteins_in_Synteny_Nb')],
                                                                line[headers_to_print.index('NetSyn_FamilyID')]))
        # COM: writting intra cluster synthesis files
        with open('{}_intraCluster_Families_NetSyn.tsv'.format(os.path.join(synthesisDirectory, cm)), 'w') as file:
            file.write('{}'.format('\t'.join(headers_to_print)))
            file.write('\n')
            for line in sorted_lines:
                file.write('{}'.format('\t'.join([str(value) for value in line])))
                file.write('\n')

    # COM: inter cluster process
    if not inter_cluster.keys():
        logger.info('This analysis does not contain edges inter clusters')
    else:
        headers_to_print = ['NetSynt_FamilyID', 'NetSyn_Clusters_Nb', 'NetSyn_Cluster_IDs', 'Target_Family', 'Syntenies_Nb', 'Species_Nb', 'Strains_Nb', 'Organisms_Nb', 'Organisms_List', 'Proteins_in_Synteny_Nb', 'Products', 'Gene_Names', 'EC_Numbers', 'Locus_Tags', 'Proteins_AC']
        if metadataFile:
            metadataHeaders = add_metadataHeaders(headers_to_print, headersMD)
            headers_to_print.extend(metadataHeaders)
        for cm, cmContent in inter_cluster.items():
            lines_to_print = []
            # COM: generation of lines containing inter cluster information
            for familyID, familyContent in cmContent.items():
                line = [''] * len(headers_to_print)
                containTarget = 'Y' if familyID in families_of_targetsIdx and set(familyContent['nodes_indexes']) & set(families_of_targetsIdx[familyID]) else 'N'
                clusters = ', '.join([str(clusterID) for clusterID in sorted(set(familyContent['clusters']))])
                nbr_clusters = len(list(set(familyContent['clusters'])))
                strains_nb = len(set(familyContent['strains']))
                species_nb = len(set(familyContent['species']))
                orgs_list = ' / '.join(familyContent['organisms'])
                orgs_nb = len(familyContent['organisms'])

                line[headers_to_print.index('NetSynt_FamilyID')] = familyID
                line[headers_to_print.index('Target_Family')] = containTarget
                line[headers_to_print.index('NetSyn_Clusters_Nb')] = nbr_clusters
                line[headers_to_print.index('NetSyn_Cluster_IDs')] = clusters
                line[headers_to_print.index('Syntenies_Nb')] = familyContent['synteny_size']
                line[headers_to_print.index('Species_Nb')] = species_nb
                line[headers_to_print.index('Strains_Nb')] = strains_nb
                line[headers_to_print.index('Organisms_Nb')] = orgs_nb
                line[headers_to_print.index('Organisms_List')] = orgs_list
                line[headers_to_print.index('Proteins_in_Synteny_Nb')] = len(set(familyContent['protein_AC'])) # use of set() no longer useful
                line[headers_to_print.index('Products')] = get_sorted_nr_data(familyContent['products'], counter=True)
                line[headers_to_print.index('Gene_Names')] = get_sorted_nr_data(familyContent['gene_names'], counter=True)
                line[headers_to_print.index('EC_Numbers')] = get_sorted_nr_data(familyContent['ec_numbers'], counter=True)
                line[headers_to_print.index('Locus_Tags')] = get_sorted_nr_data(familyContent['locus_tags'], counter=False)
                line[headers_to_print.index('Proteins_AC')] = get_sorted_nr_data(familyContent['protein_AC'], counter=False)

                if metadataFile:
                    line = add_metadata_info(line, metadataHeaders, headers_to_print, familyContent)
                lines_to_print.append(line)
            # COM: ordering lines
            order_targetFam = ['Y', 'N']
            sorted_lines = sorted(lines_to_print, key=lambda line: (-line[headers_to_print.index('NetSyn_Clusters_Nb')],
                                                                     order_targetFam.index(line[headers_to_print.index('Target_Family')]),
                                                                     -line[headers_to_print.index('Syntenies_Nb')],
                                                                     -line[headers_to_print.index('Proteins_in_Synteny_Nb')],
                                                                     sorted(line[headers_to_print.index('NetSyn_Cluster_IDs')])))
            # COM: writting inter cluster synthesis files
            with open('{}_interCluster_Families_NetSyn.tsv'.format(os.path.join(synthesisDirectory, cm)), 'w') as file:
                file.write('{}'.format('\t'.join(headers_to_print)))
                file.write('\n')
                for line in sorted_lines:
                    file.write('{}'.format('\t'.join([str(value) for value in line])))
                    file.write('\n')

    logger.info('{} completed!'.format(boxName))
    reportingMessages.append('Aller faut trouver un truc a mettre...'.format())
    common.reportingFormat(logger, boxName, reportingMessages)

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                     usage='''DataExport.py -in <inputNodes> -ie <inputEdges> -io <inputOrganisms> -ip <inputProteins> -o <outputName>\n\
\t\t[-md <metadataFile>]\n\t\t[-rrl <RedundancyRemovalLabel> | -rrt <RedundancyRemovalTaxonomy> -cm <ClusteringMethod>]''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-in', '--inputNodes', type=str,
                        required=True, help='Path of the nodes file obtained from the SyntenyFinder part')
    group1.add_argument('-ie', '--inputEdges', type=str,
                        required=True, help='Path of the edges file obtained from the SyntenyFinder part')
    group1.add_argument('-io', '--inputOrganisms', type=str,
                        required=True, help='Path of the organims file from the ClusteringIntoFamilies part')
    group1.add_argument('-ip', '--inputProteins', type=str,
                        required=True, help='Path of the proteins file from the ClusteringIntoFamilies part')
    group1.add_argument('-o', '--outputName', type=str, required=True,
                        help='Output name files')
    group1.add_argument('-md', '--metadataFile', type=str,
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
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('graphML', '{}_Results.graphML'.format(args.outputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('html', '{}_Results.html'.format(args.outputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('report', '{}_{}_report.txt'.format(args.outputName, boxName))
    #######
    # Run #
    #######
    run(args.inputNodes,
        args.inputEdges,
        args.inputOrganisms,
        args.inputProteins,
        args.metadataFile,
        args.RedundancyRemovalLabel,
        args.RedundancyRemovalTaxonomy,
        args.ClusteringMethod)
