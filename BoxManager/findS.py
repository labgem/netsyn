#!/usr/bin/env python3

##########
# Import #
##########
import argparse
import os
import logging
import pickle
import json
import sys
import math
import igraph as ig
import cairocffi
#############
# Functions #
#############

def set_userGC_similarityContext(targets, cds_info, params):
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
    for idx, afam in enumerate(target_info['similarityContext']):
        tmp_target_families.setdefault(afam, []).append(idx)
    return tmp_target_families

def intersect_families(keysA, keysB):
    return list(keysA & keysB)


def cartesian_product(listA, listB):
    cartesian_prod = [(a, b) for a in listA for b in listB]
    return cartesian_prod

def add_synton_of_targets(targetA, targetB, tmp_dict, targets_syntons):
    tarA_pos = tmp_dict[targetA]['target_pos']
    tarB_pos = tmp_dict[targetB]['target_pos']
    targets_syntons[(targetA, targetB)]['synton_of_targets'] = (tarA_pos, tarB_pos)
    if (tarA_pos, tarB_pos) in targets_syntons[(targetA, targetB)]['syntons']:
        targets_syntons[(targetA, targetB)]['boolean_targets_synton'] = True
    else:
        targets_syntons[(targetA, targetB)]['syntons'].append((tarA_pos, tarB_pos))
        targets_syntons[(targetA, targetB)]['boolean_targets_synton'] = False
    return targets_syntons

def evaluate_proximity(syntons, graph, params, look_at):
    gap = params['GAP']
    for idx, synton1 in enumerate(syntons[:-1]):
        synton2 = syntons[idx+1]
        if synton2[look_at]-synton1[look_at] <= gap+1:
            graph.add_edge(graph.vs['name'].index(synton1), graph.vs['name'].index(synton2))
    # print('details of the graph: {}'.format(graph.summary()))
    return 0

def get_connected_components(graph, synton_of_targets, params, mode='A'):
    # print('enter in get_connected_components function for the genome {}'.format(mode))
    if mode == 'A':
        look_at = 0
    else:
        look_at = 1
    sorted_syntons = sorted(graph.vs['name'], key=lambda synton: synton[look_at])
    evaluate_proximity(sorted_syntons, graph, params, look_at)
    cc = graph.components()
    for cluster in cc:
        if synton_of_targets in graph.vs[cluster]['name']:
            cc_with_target_target = cluster
            break
    # largest_cluster = []
    # for cluster in cc:
    #     if len(cluster) > len(largest_cluster):
    #         largest_cluster = cluster
    # [cc_with_target_target] = [cluster
    #                            for cluster in cc
    #                            if synton_of_targets in graph.vs[cluster]['name']]
    # print('leave get_connected_components function')
    return cc_with_target_target # largest_cluster #cc_with_target_target

def compute_score(syntons, synton_of_targets, boolean_synton_of_targets):
    minGA = min([posA for posA, posB in syntons])
    maxGA = max([posA for posA, posB in syntons])
    minGB = min([posB for posA, posB in syntons])
    maxGB = max([posB for posA, posB in syntons])

    # minGA = min([posA for posA, _ in syntons])
    # maxGA = max([posB for _, posB in syntons])
    # minGB = min([posA for posA, _ in syntons])
    # maxGB = max([posB for _, posB in syntons])
    # minGA, minGB = min([(posA[0], posB[1] for posA in syntons for posB in syntons])
    # print(minGA, maxGA, minGB, maxGB)
    synt_coverage = (maxGA-minGA+1) + (maxGB-minGB+1)
    if not boolean_synton_of_targets:
        syntons.remove(synton_of_targets)
    geneA_in_synt = len(set([synton[0] for synton in syntons]))
    geneB_in_synt = len(set([synton[1] for synton in syntons]))
    avg_genes = (geneA_in_synt + geneB_in_synt)/2
    weighting = (geneA_in_synt + geneB_in_synt)/synt_coverage
    score = avg_genes * weighting
    return score

def find_common_connected_components(maxiG, gA, gB, targetA, targetB, AB_targets_syntons, params):
    # print('enter in find_common_connected_components function')
    synton_of_targets = AB_targets_syntons['synton_of_targets']
    boolean_synton_of_targets = AB_targets_syntons['boolean_targets_synton']
    ccA = get_connected_components(gA, synton_of_targets, params)
    ccB = get_connected_components(gB, synton_of_targets, params, mode='B')
    itrsect = list(set(ccA) & set(ccB))
    # print(gB.vs[itrsect]['name'])
    # if boolean_synton_of_targets:
    if (len(itrsect) >= 2 and boolean_synton_of_targets) or (not boolean_synton_of_targets and len(itrsect) >= 3):
        # print('are long enough !')
            # print(ccA, ccB, itrsect)
        if ccA == ccB == itrsect:
            # print('we have a ccc !!!')
            score = compute_score(gA.vs[itrsect]['name'], synton_of_targets, boolean_synton_of_targets)
                # print(score)
            if not maxiG.vertex_attributes() or targetA not in maxiG.vs['name']:
                maxiG.add_vertex(targetA)
            if not maxiG.vertex_attributes() or targetB not in maxiG.vs['name']:
                maxiG.add_vertex(targetB)
            # print(targetA, targetB)
            vertex_idx_targetA = maxiG.vs['name'].index(targetA)
            vertex_idx_targetB = maxiG.vs['name'].index(targetB)
            maxiG.add_edge(vertex_idx_targetA, vertex_idx_targetB)
            edge_idx_AB = maxiG.get_eid(vertex_idx_targetA, vertex_idx_targetB)
            maxiG.es[edge_idx_AB]['weight'] = score
        else:
            # print('we need more steps')
            gA_memory = gA.copy()
            gA = ig.Graph()
            gA.add_vertices(gA_memory.vs[itrsect]['name'])
            gB = gA.copy()
            find_common_connected_components(maxiG, gA, gB, targetA, targetB, AB_targets_syntons, params)
    # else:
    #     print('not long enough :\'{')
    #     print(boolean_synton_of_targets, len(itrsect), targetA, targetB)

    # else:
    #     if len(itrsect) > 2:
    #         # print('are long enough !')
    #         # print(ccA, ccB, itrsect)
    #         if ccA == ccB == itrsect:
    #             #print('we have a ccc !!!')
    #             score = compute_score(gA.vs[itrsect]['name'], synton_of_targets, boolean_synton_of_targets)
    #             maxiG.add_vertex((targetA, targetB))
    #         else:
    #             #print('we need more steps')
    #             gA_memory = gA.copy()
    #             gA = ig.Graph()
    #             gA.add_vertices(gA_memory.vs[itrsect]['name'])
    #             gB = gA.copy()
    #             find_common_connected_components(maxiG, gA, gB, targetA, targetB, AB_targets_syntons, params)
    #     #else:
    #         #print('not long enough :\'{')
    return maxiG

def build_maxi_graph(maxiG, targets_syntons, params):
    for (targetA, targetB) in targets_syntons:
        # print('######################################')
        # print(targetA, targetB)
        gA = ig.Graph()
        gA.add_vertices(targets_syntons[(targetA, targetB)]['syntons'])
        gB = gA.copy()
        # print('first search')
        maxiG = find_common_connected_components(maxiG, gA, gB, targetA, targetB, targets_syntons[(targetA, targetB)], params)
        # if maxiG.vertex_attributes():
        #     print(maxiG.vs['name'])
    return maxiG

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

def run(gcUser, gap, gcFile, TMPDIRECTORY):
    '''
    '''
    BOXNAME = 'SyntenyFinder'
    TMPDIRECTORYPROCESS = '{}/{}'.format(TMPDIRECTORY, BOXNAME)
    if not os.path.isdir(TMPDIRECTORYPROCESS):
        os.mkdir(TMPDIRECTORYPROCESS)

    params = {
        'MAX_GC': 11,
        'USER_GC': int(gcUser),
        'GAP': int(gap),
        'INC_TARGETS_PAIR': 0
        }

    with open(gcFile, 'rb') as file:
        cds_info = pickle.load(file)
    with open('{}/new_NS_0/TMP/{}/{}'.format(TMPDIRECTORY, 'ClusteringIntoFamilies', 'targets_list.pickle'), 'rb') as file:
    #with open('{}/test/TMP/{}/{}'.format(TMPDIRECTORY, 'ClusteringIntoFamilies', 'targets_list.pickle'), 'rb') as file:
        targets_list = pickle.load(file)

    tmp_dict = {}
    targets_syntons = {}
    targets_list = sorted(list(set(targets_list))) # ligne à supprimer quand problème des doublons réglé dans CIF.py

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
                        targets_syntons.setdefault((targetA, targetB),
                                                   {}).setdefault('syntons',
                                                                  []).extend(syntons)
                    add_synton_of_targets(targetA, targetB, tmp_dict, targets_syntons)
                    # print('dictionary of targets {} and {}:\n{}\n'.format(targetA, targetB, targets_syntons[(targetA, targetB)]))

    maxi_graph = ig.Graph()
    maxi_graph = build_maxi_graph(maxi_graph, targets_syntons, params)
    # print(maxi_graph)
    # print('\n')

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
    print('\n*** Walktrap clustering ***')
    graph_walktrap = maxi_graph.community_walktrap(weights=maxi_graph.es['weight'])
    # print(graph_walktrap) # VertexDendrogram object
    walktrap_clustering = graph_walktrap.as_clustering()
    print(walktrap_clustering) # list of VertexClustering objects

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

    #layout = maxi_graph.layout('circle')
    #layout = maxi_graph.layout('drl')
    #layout = maxi_graph.layout('large')
    #layout = maxi_graph.layout('random')
    #layout = maxi_graph.layout('tree')
    #layout = maxi_graph.layout('circular')

    # layout = maxi_graph.layout('fr') # pourquoi pas ...
    #layout = maxi_graph.layout('kk') # pourquoi pas ...
    
    # pal = ig.drawing.colors.ClusterColoringPalette(len(walktrap_clustering))
    # visual_style = {
    #     'edge_width': maxi_graph.es['weight'],
    #     #'vertex_label': maxi_graph.vs['name'],
    #     'vertex_color': pal.get_many(walktrap_clustering.membership),
    #     'vertex_size': 4
    #     }
    # ig.plot(maxi_graph, "walktrap_cluster_5_3.png", layout=layout, **visual_style)


    # pal = ig.drawing.colors.ClusterColoringPalette(len(graph_louvain))
    # visual_style={
    #     'edge_width': maxi_graph.es['weight'],
    #     #'vertex_label': maxi_graph.vs['name'],
    #     'vertex_color': pal.get_many(graph_louvain.membership),
    #     'vertex_size': 4
    #     }
    # ig.plot(maxi_graph, "louvain_cluster_5_3.png", layout=layout, **visual_style)

    for cluster in range(len(walktrap_clustering)):
        for vertex in walktrap_clustering[cluster]:
            maxi_graph.vs[vertex]['cluster'] = cluster

    for cds_ref in maxi_graph.vs:
        cds_ref['protein_id'] = cds_info[cds_ref['name']]['uniprot']

    maxi_graph.write_graphml('maxi_graph.graphml')

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3], os.path.abspath('.'))
