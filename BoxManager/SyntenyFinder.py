#!/usr/bin/env python3

##########
# Import #
##########
import argparse
import os
import logging
#import pickle
#import json
#import sys
import math
import igraph as ig
#import cairocffi
import common

#############
# Functions #
#############
def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='''Description of the SyntenyFinder usage''',
                                     epilog='''All's well that ends well.''',
                                     usage='''SyntenyFinder options...''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', type=str,
                        required=True, help='Path of the input obtained from the ClusteringIntoFamilies part')
    parser.add_argument('-ws', '--WindowSize', type=int,
                        default=11, help='Window size of genomic contexts to compare (target gene inclued).\nDefault value: 11.')
    parser.add_argument('-sg', '--SyntenyGap', type=int, default=3,
                        help='Number of genes allowed betwenn tow genes in synteny.\nDefault value: 3.')
    parser.add_argument('-ssc', '--SyntenyScoreCuttoff', type=float,
                        default=0, help='Define the minimum Synteny Score Cuttoff to conserved.\nDefault value: >= 0.')
    parser.add_argument('-pn', '--ProjectName', type=str, required=True,
                        help='The project name.')
    parser.add_argument('-tl', '--TargetsList', type=str, required=True,
                        help='Path of the target list obtained from the ClusteringIntoFamilies part')
    return parser

def set_userGC_similarityContext(targets, cds_info, params):
    ''' reduce the genomic context to the user parameter '--WindowSize'
    input: list of targets, '-ws' user parameter
    output: cds_info dictionary updated with 2 new fields 1)cds_info[target]['userGC'] and
    2)cds_info[target]['similarityContext']
    '''
    half_user_window = math.floor(params['USER_GC']/2)
    for target in targets:
        center = cds_info[target]['context'].index(target)
        indices = [idx
                   for idx in range(max(0,
                                        center-half_user_window),
                                    min(len(cds_info[target]['context']),
                                        center+half_user_window+1))]
        indices = list(set(indices)) # why ???
        cds_info[target]['userGC'] = [cds_info[target]['context'][idx] for idx in indices]
        # cds_info[target]['userGC'] = [idx for idx in cds_info[target]['context'][idx] # should replace idx by cds
        #                               if idx in range(max(0,
        #                                                   center-half_user_window),
        #                                               min(len(cds_info[target]['context']),
        #                                                   center+half_user_window+1))]
        cds_info[target]['similarityContext'] = [cds_info[cds_info[target]['context'][idx]]['similarityFamily'] for idx in indices]
    return cds_info

def get_families_and_pos_in_context(target_info, tmp_target_families):
    ''' get list index of every family in a context
    input: cds_info[target] as target_info, and tmp_dict[target]['families'] as
    tmp_target_families
    output: tmp_target_families filled with a list of indices for every family
    represented in the genomic context of the target
    '''
    for idx, afam in enumerate(target_info['similarityContext']):
        tmp_target_families.setdefault(afam, []).append(idx)
    return tmp_target_families

def intersect_families(keysA, keysB):
    ''' return a list of families existing in genomic contexts of targets A AND B
    '''
    return list(keysA & keysB)


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

def add_synton_of_targets(targetA, targetB, tmp_dict, targets_syntons):
    ''' determine if the synton of targets A and B is in synteny and add it to
    the set of syntons if necessary
    input: reference of targets A and B, tmp_dict and targets_syntons dictionary
    output: updated targets_syntons where 1)'synton_of_targets' between A and B is
     characterized, 2)boolean value is set to False if this synton is not
     stored in the 'syntons' list yet, and 3)is added to the list in that case
    '''
    tarA_pos = tmp_dict[targetA]['target_pos']
    tarB_pos = tmp_dict[targetB]['target_pos']
    targets_syntons[(targetA, targetB)]['synton_of_targets'] = (tarA_pos, tarB_pos)
    if (tarA_pos, tarB_pos) in targets_syntons[(targetA, targetB)]['syntons']:
        targets_syntons[(targetA, targetB)]['boolean_targets_synton'] = True
    else:
        targets_syntons[(targetA, targetB)]['syntons'].append((tarA_pos, tarB_pos))
        targets_syntons[(targetA, targetB)]['boolean_targets_synton'] = False
    return targets_syntons

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
    return 0

def get_connected_components(graph, synton_of_targets, gapValue, mode='A'):
    ''' get the connected component of the A/B graph containing the
    'synton_of_targets'
    input: A or B graph, gapValue and the mode that allows which part of
    nodes(syntons) to look at (graph A/B -> look at index 0/1 of nodes)
    output: iGraph VertexClustering object with the 'synton_of_targets'
    '''
    if mode == 'A':
        look_at = 0
    else:
        look_at = 1
    sorted_syntons = sorted(graph.vs['name'], key=lambda synton: synton[look_at])
    evaluate_proximity(sorted_syntons, graph, gapValue, look_at)
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
    ccA = get_connected_components(gA, synton_of_targets, params['GAP'])
    ccB = get_connected_components(gB, synton_of_targets, params['GAP'], mode='B')
    itrsect = list(set(ccA) & set(ccB))
    if (len(itrsect) >= 2 and boolean_synton_of_targets) or (not boolean_synton_of_targets and len(itrsect) >= 3):
        if ccA == ccB == itrsect:
            score = compute_score(gA.vs[itrsect]['name'], synton_of_targets, boolean_synton_of_targets)
            if not maxiG.vertex_attributes() or targetA not in maxiG.vs['name']:
                maxiG.add_vertex(targetA)
                maxiG.vs[maxiG.vs['name'].index(targetA)]['targetPosition'] = synton_of_targets[0]
            if not maxiG.vertex_attributes() or targetB not in maxiG.vs['name']:
                maxiG.add_vertex(targetB)
                maxiG.vs[maxiG.vs['name'].index(targetB)]['targetPosition'] = synton_of_targets[1]
            vertex_idx_targetA = maxiG.vs['name'].index(targetA)
            vertex_idx_targetB = maxiG.vs['name'].index(targetB)
            maxiG.add_edge(vertex_idx_targetA, vertex_idx_targetB)
            edge_idx_AB = maxiG.get_eid(vertex_idx_targetA, vertex_idx_targetB)
            maxiG.es[edge_idx_AB]['weight'] = score
            params['INC_TARGETS_PAIR'] += 1
        else:
            gA_memory = gA.copy()
            gA = ig.Graph()
            gA.add_vertices(gA_memory.vs[itrsect]['name'])
            gB = gA.copy()
            find_common_connected_components(maxiG, gA, gB, targetA, targetB, AB_targets_syntons, params)
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
    logger.info('Synteny graph in construction')
    for (targetA, targetB) in targets_syntons:
        gA = ig.Graph()
        gA.add_vertices(targets_syntons[(targetA, targetB)]['syntons'])
        # COM: gA and gB have the same nodes, only edges diff
        gB = gA.copy()
        maxiG, params = find_common_connected_components(maxiG, gA, gB, targetA, targetB, targets_syntons[(targetA, targetB)], params)
    return maxiG, params

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

def run(GENOMICCONTEXTS, TARGETS_LIST, GCUSER, GAP):
    '''
    '''
    # Constants
    boxName = common.global_dict['boxName']['SyntenyFinder']
    tmpDirectoryProcess = '{}/{}'.format(common.global_dict['tmpDirectory'], boxName)
    # Outputs
    nodesOut = common.global_dict['files']['SyntenyFinder']['nodes']
    edgesOut = common.global_dict['files']['SyntenyFinder']['edges']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    logger.info('{} running...'.format(boxName))
    # Process
    if not os.path.isdir(tmpDirectoryProcess):
        os.mkdir(tmpDirectoryProcess)

    params = {
        'MAX_GC': common.global_dict['maxGCSize'],
        'USER_GC': GCUSER,
        'GAP': GAP,
        'INC_TARGETS_PAIR': 0,
        'INC_NO_SYNTENY': 0
        }

    cds_info = common.read_pickle(GENOMICCONTEXTS)
    targets_list = common.read_pickle(TARGETS_LIST)

    tmp_dict = {}
    targets_syntons = {}
    no_synteny = 0
    targets_list = sorted(list(set(targets_list))) # ligne à supprimer quand problème des doublons réglé dans CIF.py
    logger.debug('Length of the targets list: {}'.format(len(targets_list)))
    # COM: addition of last information relative to the user window size to the cds_info dictionary
    cds_info = set_userGC_similarityContext(targets_list, cds_info, params)
    for target in targets_list:
        tmp_dict[target] = {'families': {}}
        tmp_dict[target]['families'] = get_families_and_pos_in_context(cds_info[target], tmp_dict[target]['families'])
        tmp_dict[target]['target_pos'] = cds_info[target]['userGC'].index(target)

    for idx, targetA in enumerate(targets_list[:-1]):
        for targetB in targets_list[idx+1:]:
            if cds_info[targetA]['genome'] != cds_info[targetB]['genome']:
                families_intersect = intersect_families(tmp_dict[targetA]['families'].keys(),
                                                        tmp_dict[targetB]['families'].keys())
                if len(families_intersect) > 1:
                    for afam in families_intersect:
                        A_in_synt = tmp_dict[targetA]['families'][afam]
                        B_in_synt = tmp_dict[targetB]['families'][afam]
                        syntons = cartesian_product(A_in_synt, B_in_synt)
                        targets_syntons.setdefault((targetA, targetB), {}).setdefault('syntons', []).extend(syntons)
                    targets_syntons[(targetA, targetB)]['families_intersect'] = families_intersect
                    add_synton_of_targets(targetA, targetB, tmp_dict, targets_syntons)
                    # print('dictionary of targets {} and {}:\n{}\n'.format(targetA, targetB, targets_syntons[(targetA, targetB)]))
                else:
                    no_synteny += 1
            else:
                no_synteny += 1
    #logger.info('Couples of targets computed depending on synteny results')

    common.write_json(cds_info, '{}/genomicContextUser.json'.format(tmpDirectoryProcess))

    maxi_graph = ig.Graph()
    maxi_graph, params = build_maxi_graph(maxi_graph, targets_syntons, params)
    logger.info('Number of pairs of targets that don\'t share more than 1 family: {}'.format(no_synteny))
    logger.info('Number of pairs where synteny doesn\'t respect gap parameter or on target filter: {}'.format(params['INC_NO_SYNTENY']))
    logger.info('Number of pairs of targets in synteny: {}'.format(params['INC_TARGETS_PAIR']))

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
    logger.info('\n*** Walktrap clustering ***')
    graph_walktrap = maxi_graph.community_walktrap(weights=maxi_graph.es['weight'])
    walktrap_clustering = graph_walktrap.as_clustering()
    logger.info(walktrap_clustering) # list of VertexClustering objects

    for cluster in range(len(walktrap_clustering)):
        for vertex in walktrap_clustering[cluster]:
            maxi_graph.vs[vertex]['cluster_WT'] = cluster

    # ### Louvain clustering
    # print('\n*** Louvain clustering ***')
    # graph_louvain = maxi_graph.community_multilevel()
    # print(graph_louvain) # list of VertexClustering objects

    # for cluster in graph_louvain:
    #     for elt in cluster:
    #         if maxi_graph.vs[elt]['name'] not in targets_list:
    #             print(elt)
    #             print(targets_list)
    #         else:
    #             print('element {} corresponding to the target {} in the list'.format(elt, maxi_graph.vs[elt]['name']))

    list_of_nodes = []
    for target_node in maxi_graph.vs:
        cds_inc = target_node['name']
        dico = {'cds_ref': cds_inc, ### je ne sais pas encore si c'est indispensable ...
                'UniProtAC': cds_info[cds_inc]['uniprot'],
                'ProteinID': cds_info[cds_inc]['protein_id'],
                'targetPosition': maxi_graph.vs[target_node.index]['targetPosition'],
                'GC_size': len(cds_info[cds_inc]['userGC']),
                'Product': cds_info[cds_inc]['product'], ### est-ce utile, l'information sera répétée dans families
                'Contig': cds_info[cds_inc]['contig'],
                'Clustering': {'WalkTrap': maxi_graph.vs[target_node.index]['cluster_WT']},
                'Families': {}
                }
        for afam in tmp_dict[cds_inc]['families']:
            dico['Families'][afam] = {}
            dico['Families'][afam]['positions'] = tmp_dict[cds_inc]['families'][afam]
            dico['Families'][afam]['id'] = []
            dico['Families'][afam]['Protein_id'] = []
            dico['Families'][afam]['Products'] = []
            dico['Families'][afam]['EC_numbers'] = []
            for pos in dico['Families'][afam]['positions']:
                dico['Families'][afam]['id'].append(cds_info[cds_info[cds_inc]['userGC'][pos]]['uniprot'])
                dico['Families'][afam]['Protein_id'].append(cds_info[cds_info[cds_inc]['userGC'][pos]]['protein_id'])
                dico['Families'][afam]['Products'].append(cds_info[cds_info[cds_inc]['userGC'][pos]]['product'])
                dico['Families'][afam]['EC_numbers'].append(cds_info[cds_info[cds_inc]['userGC'][pos]]['ec_number'])
        list_of_nodes.append(dico)

    list_of_edges = []
    for edge in maxi_graph.es:
        targetA = min(maxi_graph.vs[edge.tuple[0]]['name'], maxi_graph.vs[edge.tuple[1]]['name'])
        targetB = max(maxi_graph.vs[edge.tuple[0]]['name'], maxi_graph.vs[edge.tuple[1]]['name'])
        if targets_syntons[targetA, targetB]:
            families = targets_syntons[(targetA, targetB)]['families_intersect']
        else:
            logger.debug('The order is not respected; the pair of targetA {} - targetB {} is referenced in the reverse order'.format(targetA, targetB))
            families = targets_syntons[(targetB,
                                        targetA)
                                       ]['families_intersect']
        dico = {'source': cds_info[targetA]['uniprot'],
                'target': cds_info[targetB]['uniprot'],
                'families': families,
                'score': maxi_graph.es[edge.index]['weight']
                }
        list_of_edges.append(dico)

    common.write_pickle(list_of_nodes, nodesOut)
    common.write_pickle(list_of_edges, edgesOut)
    common.write_json(list_of_nodes, '{}/{}'.format(tmpDirectoryProcess, 'nodes_list.json'))
    common.write_json(list_of_edges, '{}/{}'.format(tmpDirectoryProcess, 'edges_list.json'))

if __name__ == '__main__':
    parser = argumentsParser()
    args = parser.parse_args()
    BOXNAME = 'SyntenyFinder'
    TMPDIRECTORY = '{}/TMP'.format(args.ProjectName)
    MAXGCSIZE = 11
    run(BOXNAME, TMPDIRECTORY, args.input, args.TargetsList, MAXGCSIZE, args.WindowSize, args.SyntenyGap, args.SyntenyScoreCuttoff)
