#!/usr/bin/env python3

##########
# Import #
##########
import argparse
import os
import logging
import math
import igraph as ig
import shutil
import markov_clustering as mc
import networkx as nx
import common

#############
# Functions #
#############

# def get_userGC(targets_info, windowSize, prots_info):
#     ''' reduce the genomic context to the user parameter '--WindowSize'
#     input: list of targets_info, '-ws' user parameter
#     output: prots_info dictionary updated with 2 new fields 1)prots_info[target]['userGC'] and
#     2)prots_info[target]['similarityContext']
#     '''
#     half_user_window = math.floor(windowSize/2)
#     new_targets_info = {}

#     proteinsToConserve = []
#     proteinsIndexToDel = []
#     for target_idx, target_dict in targets_info.items():
#         center = target_dict['context_idx'].index(target_idx)
#         low_limit = max(0, center-half_user_window)
#         high_limit = min(len(target_dict['context_idx']),
#                          center+half_user_window+1)

#         for i in target_dict['context_idx'][:low_limit] + target_dict['context_idx'][high_limit:]:
#             if i not in proteinsToConserve and i not in proteinsIndexToDel:
#                 proteinsIndexToDel.append(i)
#         for i in target_dict['context_idx'][low_limit:high_limit]:
#             if i in proteinsIndexToDel:
#                 del proteinsIndexToDel[proteinsIndexToDel.index(i)]
#             if i not in proteinsToConserve:
#                 proteinsToConserve.append(i)
#     proteinsIndexes = [i for i in range(len(prots_info))]

#     decrement = 0
#     for i in proteinsIndexes:
#         if i in proteinsIndexToDel:
#             proteinsIndexes[i] = -1
#             decrement += 1
#         else:
#             proteinsIndexes[i] -= decrement

#     for target_idx, target_dict in targets_info.items():#
#         center = target_dict['context_idx'].index(target_idx)#
#         low_limit = max(0, center-half_user_window)#
#         high_limit = min(len(target_dict['context_idx']),
#                          center+half_user_window+1)#
#         target_dict['context_idx'] = [proteinsIndexes[i] for i in target_dict['context_idx'][low_limit:high_limit]]# juste i
#         target_dict['context'] = target_dict['context'][low_limit:high_limit]#
#         target_dict['target_pos'] = target_dict['context'].index(target_dict['id'])#
#         new_index = target_dict['context_idx'][target_dict['target_pos']]
#         new_targets_info.setdefault(new_index, {}).update(target_dict)

#     for i in sorted(proteinsIndexToDel, reverse=True):
#         del prots_info[i]
#     return new_targets_info, prots_info

def get_userGC(targets_info, windowSize):
    ''' reduce the genomic context to the user parameter '--WindowSize'
    input: list of targets_info, '-ws' user parameter
    output: prots_info dictionary updated on 2 fields 1)prots_info[target]['context_idx'] and
    2)prots_info[target]['context']
    '''
    half_user_window = math.floor(windowSize/2)
    for target_idx, target_dict in enumerate(targets_info):
        center = target_dict['context_idx'].index(target_dict['protein_idx'])
        low_limit = max(0, center-half_user_window)
        high_limit = min(len(target_dict['context_idx']),
                         center+half_user_window+1)
        targets_info[target_idx]['context_idx'] = [i for i in target_dict['context_idx'][low_limit:high_limit]]
        targets_info[target_idx]['context'] = target_dict['context'][low_limit:high_limit]
        targets_info[target_idx]['target_pos'] = target_dict['context'].index(target_dict['id'])
    return targets_info

def proteinsRemoval(prots_info, targets_info, maxi_graph):#, targets_syntons):
    '''
    Removes:
        I. the proteins exclued by le resizing of genomic context
       II. the proteins of genomic contexts for the targets without the conserved synteny.

    Updates the protein indexes in prots_info, targets_info, maxi_graph, targets_syntons.
    '''
    # Defines the proteins to delete
    proteinsIndex = [i for i in range(len(prots_info))]
    proteinsIndexInSynteny = []
    #print(len(maxi_graph.vs['name']))
    for targetIndex in maxi_graph.vs['name']:
        for proteinIndex in targets_info[targetIndex]['context_idx']:
            if proteinIndex not in proteinsIndexInSynteny:
                proteinsIndexInSynteny.append(int(proteinIndex))
    proteinsIndexInSynteny.sort()
    proteinsIndexToDel = sorted(list(set(proteinsIndex).difference(set(proteinsIndexInSynteny))))

    # Defines the new indexes. (the poroteins indexes of proteins to delete becomes -1, the others is decremente.)
    decrement = 0
    for i in proteinsIndex:
        if i in proteinsIndexToDel:
            proteinsIndex[i] = -1
            decrement += 1
        else:
            proteinsIndex[i] -= decrement

    # Uptades indexes in targets_info
    new_targets_info = []
    for target_idx, target_dict in enumerate(targets_info):
        target_dict['context_idx'] = [str(proteinsIndex[int(i)]) for i in target_dict['context_idx']]
        new_protein_index = str(target_dict['context_idx'][target_dict['target_pos']])
        target_dict['protein_idx'] = new_protein_index
        new_targets_info.append(target_dict)
    #print(len(new_targets_info))

    # Uptades indexes in maxi_graph
    # maxi_graph.vs['name'] = [str(proteinsIndex[int(lastIndex)]) for lastIndex in maxi_graph.vs['name']]

    # Uptades indexes in new_targets_syntons
    # new_targets_syntons = {}
    # for (indexA, indexB), value in targets_syntons.items():
    #     new_targets_syntons[(str(proteinsIndex[int(indexA)]), str(proteinsIndex[int(indexB)]))] = value

    # Removes proteins that are not preserved.
    for i in sorted(proteinsIndexToDel, reverse=True):
        del prots_info[i]

    # Uptades indexes in prots_info[idx]['targets_idx']
    for protein_idx, _ in enumerate(prots_info):
        prots_info[protein_idx]['targets_idx'] = []
    for target_idx, target_dict in enumerate(new_targets_info):
        for idx in target_dict['context_idx']:
            prots_info[int(idx)]['targets_idx'].append(target_dict['protein_idx'])

    return new_targets_info, prots_info#, maxi_graph#, new_targets_syntons

def get_families_and_pos_in_context(target_dict):
    ''' get list index of every family in a context
    input: prots_info[target] as target_info, and tmp_dict[target]['families'] as
    tmp_target_families
    output: tmp_target_families filled with a list of indices for every family
    represented in the genomic context of the target
    '''
    for idx, afam in enumerate(target_dict['families']):
        target_dict.setdefault('families_positions', {}).setdefault(afam, []).append(idx)
    return target_dict

def intersect_families(keysA, keysB):
    ''' return a list of families existing in genomic contexts of targets A AND B
    '''
    return list(set(keysA) & set(keysB))


def cartesian_product(listA, listB):
    '''
    input: lists of indices stored in tmp_dict[targetA(B)]['families'][afam]
    output: list containing all possible pairs between listA and listB
    i.e.
    listA = [1, 2]
    listB = [0, 2, 6]
    output = [(1,0), (1,2), (1,6), (2,0), (2,2), (2,6)]
    '''
    cartesian_prod = [(a, b) for a in listA for b in listB]
    return cartesian_prod

def add_synton_of_targets(targetA, targetB, target_synton_AB):
    ''' determine if the synton of targets A and B is in synteny and add it to
    the set of syntons if necessary
    input: reference of targets A and B, tmp_dict and targets_syntons dictionary
    output: updated targets_syntons where 1)'synton_of_targets' between A and B is
     characterized, 2)boolean value is set to False if this synton is not
     stored in the 'syntons' list yet, and 3)is added to the list in that case
    '''
    tarA_pos = targetA['target_pos']
    tarB_pos = targetB['target_pos']
    target_synton_AB['synton_of_targets'] = (tarA_pos, tarB_pos)
    if (tarA_pos, tarB_pos) in target_synton_AB['syntons']:
        target_synton_AB['boolean_targets_synton'] = True
    else:
        target_synton_AB['syntons'].append((tarA_pos, tarB_pos))
        target_synton_AB['boolean_targets_synton'] = False
    return target_synton_AB

def evaluate_proximity(syntons, graph, gapValue, look_at):
    ''' add new edges in graph when nodes are enough close to each other
    input: list of syntons between targets A and B, specific graph of genomeA or
    genomeB, user parameter for '--SyntenyGap', and value look_at (0/1 for the
    genomic context A/B)
    output: no output, but the graph has a new edge if the gapValue constraint is fulfilled
    '''
    for idx, synton1 in enumerate(syntons[:-1]):
        synton2 = syntons[idx+1]
        if synton2[look_at]-synton1[look_at] <= gapValue+1:
            graph.add_edge(graph.vs['name'].index(synton1), graph.vs['name'].index(synton2))
    return graph

def get_connected_components(graph, synton_of_targets, gapValue, look_at):
    ''' get the connected component of the A/B graph containing the
    'synton_of_targets'
    input: A or B graph, gapValue and the mode that allows which part of
    nodes(syntons) to look at (graph A/B -> look at index 0/1 of nodes)
    output: iGraph VertexClustering object with the 'synton_of_targets'
    '''
    sorted_syntons = sorted(graph.vs['name'], key=lambda synton: synton[look_at])
    graph = evaluate_proximity(sorted_syntons, graph, gapValue, look_at)
    ccs = graph.components()
    for component_connexe in ccs:
        if synton_of_targets in graph.vs[component_connexe]['name']:
            cc_with_target_target = component_connexe
            break
    return cc_with_target_target

def compute_score(syntons, synton_of_targets, boolean_synton_of_targets):
    ''' return the score of the synteny between targets A and B
    if boolean value is false, the node is removed from connected component list
    (syntons)
    '''
    targetPosA = synton_of_targets[0]
    targetPosB = synton_of_targets[1]
    if not boolean_synton_of_targets:
        syntons.remove(synton_of_targets)
    minGA = min(targetPosA+1, min([posA for posA, posB in syntons]))
    maxGA = max(targetPosA-1, max([posA for posA, posB in syntons]))
    minGB = min(targetPosB+1, min([posB for posA, posB in syntons]))
    maxGB = max(targetPosB-1, max([posB for posA, posB in syntons]))
    ### COM: suppression of the target position if it is inbetween the extremum and not in synteny
    ### to not create a gap ponderation on the target
    if minGA < targetPosA < maxGA and targetPosA not in [synton[0] for synton in syntons]:
        syntons = [(posA-1, posB) if posA > targetPosA else (posA, posB) for posA, posB in syntons]
        maxGA = maxGA-1
    if minGB < targetPosB < maxGB and targetPosB not in [synton[1] for synton in syntons]:
        syntons = [(posA, posB-1) if posB > targetPosB else (posA, posB) for posA, posB in syntons]
        maxGB = maxGB-1
    ### COM: score calculation
    synt_coverage = (maxGA-minGA+1) + (maxGB-minGB+1)
    geneA_in_synt = len(set([synton[0] for synton in syntons]))
    geneB_in_synt = len(set([synton[1] for synton in syntons]))
    avg_genes = (geneA_in_synt + geneB_in_synt)/2
    weighting = (geneA_in_synt + geneB_in_synt)/synt_coverage
    score = avg_genes * weighting
    return score

def find_common_connected_components(maxiG, gA, gB, targetA, targetB, AB_targets_syntons, params):
    ''' recursive function to detect a common connected component set in the set
    of nodes in genome A and genome B
    input: maxi graph, subgraphs A and B (nodes are identical, edges differ),
    targets_sytons[(targetA, targetB)] dictionary as AB_targets_syntons
    output: updated params dictionary, and updated maxi graph regarding
    1)added nodes if necessary, 2) added edges weighted by the score between nodes that represent
    targets A and B
    '''
    synton_of_targets = AB_targets_syntons['synton_of_targets']
    boolean_synton_of_targets = AB_targets_syntons['boolean_targets_synton']
    ccA = get_connected_components(gA, synton_of_targets, params['GAP'], look_at=0)
    ccB = get_connected_components(gB, synton_of_targets, params['GAP'], look_at=1)
    itrsect = list(set(ccA) & set(ccB))
    if (len(itrsect) >= 2 and boolean_synton_of_targets) or (not boolean_synton_of_targets and len(itrsect) >= 3):
        if ccA == ccB == itrsect:
            score = compute_score(gA.vs[itrsect]['name'], synton_of_targets, boolean_synton_of_targets)
            if score < params['SCORE_CUTOFF']:
                params['INC_CUTOFF'] += 1
            else:
                if not maxiG.vertex_attributes() or targetA not in maxiG.vs['name']:
                    maxiG.add_vertex(targetA)
                    maxiG.vs[maxiG.vs['name'].index(targetA)]['targetPosition'] = synton_of_targets[0]
                if not maxiG.vertex_attributes() or targetB not in maxiG.vs['name']:
                    maxiG.add_vertex(targetB)
                    maxiG.vs[maxiG.vs['name'].index(targetB)]['targetPosition'] = synton_of_targets[1]
                # print(maxiG.vs['name'])
                # print("\t"*2,targetA, targetB)
                vertex_idx_targetA = maxiG.vs['name'].index(targetA)
                vertex_idx_targetB = maxiG.vs['name'].index(targetB)
                # print("\t"*2,vertex_idx_targetA, vertex_idx_targetB)
                maxiG.add_edge(vertex_idx_targetA, vertex_idx_targetB)
                edge_idx_AB = maxiG.get_eid(vertex_idx_targetA, vertex_idx_targetB)
                maxiG.es[edge_idx_AB]['weight'] = score
        else:
            gA_memory = gA.copy()
            gA = ig.Graph()
            gA.add_vertices(gA_memory.vs[itrsect]['name'])
            gB = gA.copy()
            maxiG, params = find_common_connected_components(maxiG, gA, gB, targetA, targetB, AB_targets_syntons, params)
    else:
        params['INC_NO_SYNTENY'] += 1
    return maxiG, params

def build_maxi_graph(maxiG, targets_syntons, params):
    ''' construction of the maxi graph where nodes are equivalent to targets and
    edges represent a synteny relation between targets
    input: empty graph (maxiG), targets_syntons dictionary
    output: graph with every synteny relation weighted by the score of the
    relation
    '''
    logger = logging.getLogger('{}.{}'.format(build_maxi_graph.__module__, build_maxi_graph.__name__))
    logger.info('Synteny graph in construction...')
    for (targetA, targetB) in targets_syntons:
        gA = ig.Graph()
        gA.add_vertices(targets_syntons[(targetA, targetB)]['syntons'])
        # COM: gA and gB have the same nodes, only edges diff
        gB = gA.copy()
        maxiG, params = find_common_connected_components(maxiG, gA, gB, targetA, targetB, targets_syntons[(targetA, targetB)], params)
        # print(targetA, targetB)
    return maxiG, params

def get_proteins_in_synteny(families, targetAinfo, targetBinfo, prots_info):
    '''
    *******************
    on peut faire aussi la recherche via les syntons
    pas besoin de prots_info dans les paramètres de la fonction ...
    *******************
    '''
    proteins_idx_source = []
    proteins_idx_target = []
    for afam in families:
        for position in targetAinfo['families_positions'][afam]:
            proteins_idx_source.append(targetAinfo['context_idx'][position])
        for position in targetBinfo['families_positions'][afam]:
            proteins_idx_target.append(targetBinfo['context_idx'][position])
    return proteins_idx_source, proteins_idx_target

# def fix_dendrogram(graph, cl):
#     '''known bug with incomplete dendrograms
#     https://lists.nongnu.org/archive/html/igraph-help/2014-02/msg00067.html
#     takes a graph and an incomplete dendrogram and completes the dendrogram by merging
#     the remaining nodes in arbitrary order
#     '''
#     already_merged = set()
#     for merge in cl.merges:
#         already_merged.update(merge)
#     num_dendrogram_nodes = graph.vcount() + len(cl.merges)
#     print(num_dendrogram_nodes)
#     not_merged_yet = sorted(set(range(num_dendrogram_nodes)) - already_merged)
#     print(not_merged_yet)
#     if len(not_merged_yet) < 2:
#         return
#     v1, v2 = not_merged_yet[:2]
#     cl._merges.append((v1, v2))
#     del not_merged_yet[:2]
#     missing_nodes = range(num_dendrogram_nodes,
#                           num_dendrogram_nodes + len(not_merged_yet))
#     print(missing_nodes)
#     cl._merges.extend(zip(not_merged_yet, missing_nodes))
#     cl._nmerges = graph.vcount()-1
#     return cl

def run(PROTEINS, TARGETS, GCUSER, GAP, CUTOFF, ADVANCEDSETTINGSFILENAME):
    '''
    '''
    # Constants
    boxName = common.global_dict['boxName']['SyntenyFinder']
    dataDirectoryProcess = os.path.join(common.global_dict['dataDirectory'], boxName)
    # Outputs
    nodesOut = common.global_dict['files']['SyntenyFinder']['nodes']
    edgesOut = common.global_dict['files']['SyntenyFinder']['edges']
    protsOut = common.global_dict['files']['SyntenyFinder']['proteins']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    reportingMessages = []
    print('')
    logger.info('{} running...'.format(boxName))
    # Process
    if not os.path.isdir(dataDirectoryProcess):
        os.mkdir(dataDirectoryProcess)

    params = {
        'MAX_GC': common.global_dict['maxGCSize'],
        'USER_GC': GCUSER,
        'GAP': GAP,
        'SCORE_CUTOFF': CUTOFF,
        'INC_CUTOFF': 0,
        'INC_NO_SYNTENY': 0
        }

    advanced_settings = common.readYamlAdvancedSettingsFile(ADVANCEDSETTINGSFILENAME, common.getClusteringMethodsDefaultSettings())

    prots_info = common.readJSON(PROTEINS, common.getProteinsFamiliesStepSchema())
    targets_info = common.readJSON(TARGETS, common.getTargetsTaxonomyStepschema())

    targets_syntons = {}
    no_synteny = 0
    # logger.info('Length of the targets list: {}'.format(len(targets_info)))
    # COM: addition of last information relative to the user window size to the prots_info dictionary
    if params['MAX_GC'] != params['USER_GC']:
        logger.info('Genomic context resizing ({} -> {})...'.format(params['MAX_GC'], params['USER_GC']))
        targets_info = get_userGC(targets_info, params['USER_GC'])

    for target_dict in targets_info:
        target_dict['families'] = [prots_info[int(idx)]['family'] for idx in target_dict['context_idx']]
        target_dict = get_families_and_pos_in_context(target_dict)
        if params['MAX_GC'] == params['USER_GC']:
            target_dict['target_pos'] = target_dict['context_idx'].index(target_dict['protein_idx'])

    # targets_list = list(targets_info.keys())
    # for idx, targetAidx in enumerate(targets_list[:-1]):
        # targetA = targets_info[targetAidx]
        # for targetBidx in targets_list[idx+1:]:
        #     targetB = targets_info[targetBidx]
    for targetAidx, targetA in enumerate(targets_info[:-1]):
        for i, targetB in enumerate(targets_info[targetAidx+1:]):
            targetBidx = (targetAidx + i + 1)
            if targetA['organism_idx'] != targetB['organism_idx']:
                families_intersect = intersect_families(targetA['families'],
                                                        targetB['families'])
                fam_intersect_len = len(families_intersect) if not None in families_intersect else len(families_intersect)-1
                if fam_intersect_len > 1:
                    for afam in [family for family in families_intersect if family != None]:
                        A_in_synt = targetA['families_positions'][afam]
                        B_in_synt = targetB['families_positions'][afam]
                        syntons = cartesian_product(A_in_synt, B_in_synt)
                        targets_syntons.setdefault((targetAidx, targetBidx), {}).setdefault('syntons', []).extend(syntons)
                    targets_syntons[(targetAidx, targetBidx)]['families_intersect'] = [family for family in families_intersect if family != None]
                    targets_syntons[(targetAidx, targetBidx)] = add_synton_of_targets(targetA, targetB, targets_syntons[(targetAidx, targetBidx)])
                else:
                    no_synteny += 1
            else:
                no_synteny += 1
    maxi_graph = ig.Graph()
    maxi_graph, params = build_maxi_graph(maxi_graph, targets_syntons, params)
    # logger.info('Number of pairs of targets that don\'t share more than 1 family: {}'.format(no_synteny))
    # logger.info('Number of conserved synteny between 2 targets doesn\'t respect gap parameter: {}'.format(params['INC_NO_SYNTENY']))
    # logger.info('Number of conserved synteny between 2 targets where synteny score is less than Synteny Score Cut-Off: {}'.format(params['INC_CUTOFF']))
    # logger.info('Number of conserved synteny between 2 targets in the analysis: {}'.format(len(maxi_graph.es)))

    ### Edge-betweenness clustering
    # graph_edge_btwness = maxi_graph.community_edge_betweenness(directed=False)
    # print(graph_edge_btwness)
    # complete_dendogram = fix_dendrogram(maxi_graph, graph_edge_btwness)
    # btwness_clusters = complete_dendogram.as_clustering()
    # print(btwness_clusters)

    # with open('maxi_graph_11_1.list', 'w') as file:
    #     ordered_vertices = sorted(maxi_graph.vs['name'])
    #     for vertex in ordered_vertices:
    #         file.write('{}\t{}\n'.format(maxi_graph.vs['name'].index(vertex), vertex))

    ### Walktrap clustering
    graph_walktrap = maxi_graph.community_walktrap(weights='weight', steps=advanced_settings[common.global_dict['WalkTrap']]['walktrap_step'])
    walktrap_clustering = graph_walktrap.as_clustering()

    for cluster in range(len(walktrap_clustering)):
        for vertex in walktrap_clustering[cluster]:
            maxi_graph.vs[vertex]['cluster_WT'] = cluster

    ### Louvain clustering
    graph_louvain = maxi_graph.community_multilevel(weights='weight')

    for cluster in range(len(graph_louvain)):
        for vertex in graph_louvain[cluster]:
            maxi_graph.vs[vertex]['cluster_Louvain'] = cluster

    # ### Spinglass clustering
    # # algorithm dont le résultat peut varier (changements mineurs)
    # # méthode trop lente ~40sec pour clusteriser 526 noeuds
    # logger.info('*** Spinglass clustering ***')
    # connected_components = maxi_graph.components()
    # for cluster in connected_components:
    #     graph_cluster = maxi_graph.subgraph(cluster)
    #     graph_spinglass = graph_cluster.community_spinglass()#weigths=maxi_graph.es['weight'])
    #     logger.info('{} -- {}'.format(cluster, graph_spinglass))

    # ### Label Propagation Clustering
    # # solutions variant entre les résultats retournés par Louvain et WalkTrap
    # # cependant, ne prends pas en compte le poids des arêtes
    # # algorithme stochastique, il faudrait obtenir un clustering consensus à partir de plusieurs runs (>100-200)
    # logger.info('*** Label Propagation Clustering***')
    # graph_labelPropagation = maxi_graph.community_label_propagation()#weigths=maxi_graph.es['weight'])
    # logger.info(graph_labelPropagation) # list of VertexClustering objects

    ### Infomap Clustering
    # donne le même résultat que Louvain (sur données UniProtAC, pas avec BKACE !)
    # l'algo n'utilise pas non plus le poids des arêtes
    graph_infomap = maxi_graph.community_infomap(edge_weights='weight', trials=advanced_settings[common.global_dict['Infomap']]['infomap_trials'])

    for cluster in range(len(graph_infomap)):
        for vertex in graph_infomap[cluster]:
            maxi_graph.vs[vertex]['cluster_Infomap'] = cluster

    # ### Leading EigenVector Clustering
    # # donne le même résultat que WalkTrap (sur données UniProtAC, pas avec BKACE !)
    # # ne génère que 11 clusters à partir d'un graph formé de 10 composantes connexes (en split 1 seul en 2)
    # logger.info('*** Leading EigenVector Clustering***')
    # graph_eigenvector = maxi_graph.community_leading_eigenvector()#weigths=maxi_graph.es['weight'])
    # logger.info(graph_eigenvector) # list of VertexClustering objects

    ### Markov Cluster Algorithm (MCL Clustering)
    nxGraph = nx.Graph()
    names = maxi_graph.vs['name']
    nxGraph.add_nodes_from(names)
    nxGraph.add_weighted_edges_from([(names[x[0]], names[x[1]], maxi_graph.es[maxi_graph.get_eid(x[0],x[1])]['weight']) for x in maxi_graph.get_edgelist()])

    matrix_adjacency = nx.to_scipy_sparse_matrix(nxGraph, weight='weight')

    # # A tester sur d'autres jeux de données pour voir l'évolution des paramètres inflation et expansion
    # for inflation in [i/10 for i in range(15,26)]:
    #     for expansion in [j for j in range(2,11)]:
    #         result = mc.run_mcl(matrix_adjacency, inflation=inflation, expansion=expansion)
    #         clusters = mc.get_clusters(result)
    #         Q = mc.modularity(matrix=result, clusters=clusters)
    #         print('inflation: {}\texpansion: {}\t modularity: {}'.format(inflation, expansion, Q))

    result = mc.run_mcl(
        matrix_adjacency,
        inflation=advanced_settings[common.global_dict['MCL']]['MCL_inflation'],
        expansion=advanced_settings[common.global_dict['MCL']]['MCL_expansion'],
        iterations=advanced_settings[common.global_dict['MCL']]['MCL_iterations']
    )
    clusters = mc.get_clusters(result)

    for cluster in range(len(clusters)):
        for vertex in clusters[cluster]:
            maxi_graph.vs[vertex]['cluster_MCL'] = cluster

    targetsNumber = len(targets_info)
    #print(len(targets_info), len(maxi_graph.vs['name']))
    # targets_info, prots_info, maxi_graph, targets_syntons = proteinsRemoval(prots_info, targets_info, maxi_graph, targets_syntons)
    targets_info, prots_info = proteinsRemoval(prots_info, targets_info, maxi_graph)
    #print(len(targets_info), len(maxi_graph.vs['name']))
    #print(set(targets_info.keys()).difference(set(maxi_graph.vs['name'])))
    # print(targets_info.keys())
    # print(maxi_graph.vs['name'])

    list_of_nodes = []
    for target_node in maxi_graph.vs:
        target_idx = int(target_node['name'])
        protein_idx = int(targets_info[target_idx]['protein_idx'])
        dico = {'protein_idx': protein_idx,
                'id': prots_info[protein_idx]['id'],
                'UniProt_AC': prots_info[protein_idx]['UniProt_AC'],
                'protein_AC': prots_info[protein_idx]['protein_AC'],
                #'targetPosition': maxi_graph.vs[target_node.index]['targetPosition'],
                #'GC_size': len(prots_info[cds_inc]['userGC']),
                #'Product': prots_info[cds_inc]['product'], ### est-ce utile, l'information sera répétée dans families
                'context': targets_info[target_idx]['context'],
                'context_idx': targets_info[target_idx]['context_idx'],
                'organism_id': targets_info[target_idx]['organism_id'],
                'organism_idx': targets_info[target_idx]['organism_idx'],
                'clusterings': {'WalkTrap':
                                    maxi_graph.vs[target_node.index]['cluster_WT'],
                                'Louvain':
                                    maxi_graph.vs[target_node.index]['cluster_Louvain'],
                                'Infomap':
                                    maxi_graph.vs[target_node.index]['cluster_Infomap'],
                                'MCL':
                                    maxi_graph.vs[target_node.index]['cluster_MCL']
                               },
                'families': list(set(targets_info[target_idx]['families'])),
                'Size': 1
                }
        list_of_nodes.append(dico)

    list_of_edges = []
    for edge in maxi_graph.es:
        node_0_idx = maxi_graph.vs[edge.tuple[0]].index
        node_1_idx = maxi_graph.vs[edge.tuple[1]].index
        targetA = min(int(maxi_graph.vs[node_0_idx]['name']), int(maxi_graph.vs[node_1_idx]['name']))
        targetB = max(int(maxi_graph.vs[node_0_idx]['name']), int(maxi_graph.vs[node_1_idx]['name']))
        nodeA_idx, nodeB_idx = (node_0_idx, node_1_idx) if maxi_graph.vs[node_0_idx]['name'] == targetA else (node_1_idx, node_0_idx)

        if targets_syntons[(targetA, targetB)]:
            families = targets_syntons[(targetA, targetB)]['families_intersect']
        else:
            logger.debug('The order is not respected; the pair of targetA {} - targetB {} is referenced in the reverse order'.format(targetA, targetB))
            families = targets_syntons[(targetB, targetA)]['families_intersect']

        proteins_idx_source, proteins_idx_target = get_proteins_in_synteny(families, targets_info[targetA], targets_info[targetB], prots_info)

        dico = {'source': nodeA_idx,
                'target': nodeB_idx,
                'proteins_idx_source': proteins_idx_source,
                'proteins_idx_target': proteins_idx_target,
                'weight': maxi_graph.es[edge.index]['weight']
                }
        list_of_edges.append(dico)

    common.write_json(prots_info, protsOut)
    common.write_json(list_of_nodes, nodesOut)
    common.write_json(list_of_edges, edgesOut)

    if os.listdir(dataDirectoryProcess) == []:
        shutil.rmtree(dataDirectoryProcess)

    logger.info('{} completed!'.format(boxName))
    # logger.info('Number of pairs of targets that don\'t share more than 1 family: {}'.format(no_synteny))
    # logger.info('Number of conserved synteny between 2 targets doesn\'t respect gap parameter: {}'.format(params['INC_NO_SYNTENY']))
    # logger.info('Number of conserved synteny between 2 targets where synteny score is less than Synteny Score Cut-Off: {}'.format(params['INC_CUTOFF']))
    # logger.info('Number of conserved synteny between 2 targets in the analysis: {}'.format(len(maxi_graph.es)))
    reportingMessages.append('Genomic context size: {}'.format(GCUSER))
    reportingMessages.append('Protein targets number with a conserved genomic context: {}/{}'.format(
        len(list_of_nodes),targetsNumber
    ))
    reportingMessages.append('Conserved genomic context found: {}'.format(len(list_of_edges)))
    reportingMessages.append('NetSyn Walktrap: {}'.format(walktrap_clustering.summary()))
    reportingMessages.append('NetSyn Louvain: {}'.format(graph_louvain.summary()))
    reportingMessages.append('NetSyn Infomap: {}'.format(graph_infomap.summary()))
    reportingMessages.append('NetSyn MCL: Clustering with {} elements and {} clusters'.format(len([item for subCluster in clusters for item in subCluster]), len(clusters)))
    common.reportingFormat(logger, boxName, reportingMessages)

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                     usage='''SyntenyFinder.py -ip <proteinsFileName> -it <targetsFileName> -o <outputName>\n\
\t\t-ws <WindowSize> -sg <SyntenyGap> -ssc <SyntenyScoreCutoff>''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-ip', '--inputProteins', type=str,
                        required=True, help='Proteins File')
    group1.add_argument('-it', '--inputTargets', type=str,
                        required=True, help='Targets File')
    group1.add_argument('-o', '--outputName', type=str,
                        required=True, help='Output name files')
    group1.add_argument('-ws', '--WindowSize', type=int,
                        default=common.global_dict['maxGCSize'],
                        choices=common.widowsSizePossibilities(common.global_dict['minGCSize'],common.global_dict['maxGCSize']),
                        help='Window size of genomic contexts to compare (target gene inclued).\nDefault value: {}.'.format(common.global_dict['maxGCSize']))
    group1.add_argument('-sg', '--SyntenyGap', type=int, default=3,
                        help='Number of genes allowed between two genes in synteny.\nDefault value: 3')
    group1.add_argument('-ssc', '--SyntenyScoreCutoff', type=float,
                        default=common.global_dict['sscDefault'], help='Define the minimum Synteny Score Cut off to conserved.\nDefault value: >= {}'.format(common.global_dict['sscDefault']))

    group3 = parser.add_argument_group('Advanced settings')
    group3.add_argument('-asc', '--ClusteringAdvancedSettings', type=str,
                        help='YAML file with the advanced clustering settings to determine synteny NetSyn. Settings of clusterings methods')


    group2 = parser.add_argument_group('logger')
    group2.add_argument( '--log_level',
                         type = str,
                         nargs = '?',
                         default = 'INFO',
                         help = 'log level',
                         choices = ['ERROR', 'error', 'WARNING', 'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                         required = False )
    group2.add_argument( '--log_file',
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
    common.global_dict['dataDirectory'] = '.'
    boxName = common.global_dict['boxName']['SyntenyFinder']
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault('nodes', '{}_nodes.json'.format(args.outputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault('edges', '{}_edges.json'.format(args.outputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault('proteins', '{}_proteins_syntenyStep.json'.format(args.outputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault('report', '{}_{}_report.txt'.format(args.outputName, boxName))
    #######
    # Run #
    #######
    run(args.inputProteins, args.inputTargets ,args.WindowSize, args.SyntenyGap, args.SyntenyScoreCutoff, args.ClusteringAdvancedSettings)
