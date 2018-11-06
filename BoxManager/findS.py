#!/usr/bin/env python3

##########
# Import #
##########
import argparse
import os
import logging
import pickle
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
        center = cds_info[target]['context'].index((target, cds_info[target]['protein_id']))
        indices = [idx
                   for idx in range(max(0,
                                        center-half_user_window),
                                    min(len(cds_info[target]['context']),
                                        center+half_user_window+1))]
        indices = list(set(indices)) # why ???
        cds_info[target]['userGC'] = [cds_info[target]['context'][idx][0] for idx in indices]
        cds_info[target]['similarityContext'] = [cds_info[cds_info[target]['context'][idx][0]]['similarityFamily'] for idx in indices]
    return cds_info

def intersect_families(targetA, targetB):
    familiesA = targetA['similarityContext']
    familiesB = targetB['similarityContext']
    intersect = list(set(familiesA) & set(familiesB))
    return intersect

def get_pos_synton(target, family, cds_info):
    pos_in_synt = [cds_info[target]['userGC'].index(cds)
                   for cds in cds_info[target]['userGC']
                   if cds_info[cds]['similarityFamily'] == family]
    return pos_in_synt

def cartesian_product(listA, listB):
    cartesian_prod = [(a, b) for a in listA for b in listB]
    return cartesian_prod

def get_target_pos(context, target):
    '''if comprehension list return 0 or more than 1 element,
    it raise "ValueError: too many values to unpack (expected 1)"
    useful when only one value expected (in Clustering_Into_Families)
    '''
    return [context.index((cds_ref, cds_id))
            for (cds_ref, cds_id) in context
            if cds_ref == target]

def add_target_target_synton(targetA, targetB, cds_info, targets_syntons):
    [tarA_pos] = get_target_pos(cds_info[targetA]['context'], targetA)
    [tarB_pos] = get_target_pos(cds_info[targetB]['context'], targetB)

    if (tarA_pos, tarB_pos) in targets_syntons[(targetA, targetB)]['syntons']:
        targets_syntons[(targetA, targetB)]['AB_in_synt'] = (tarA_pos, tarB_pos, True)
    else:
        targets_syntons[(targetA, targetB)]['syntons'].append((tarA_pos, tarB_pos))
        targets_syntons[(targetA, targetB)]['AB_in_synt'] = (tarA_pos, tarB_pos, False)
    return targets_syntons

def evalute_proximity(syntons, graph, params, look_at):
    gap = params['GAP']
    for idx, synton1 in enumerate(syntons[:-1]):
        synton2 = syntons[idx+1]
        if synton2[look_at]-synton1[look_at] <= gap+1:
            graph.add_edge(graph.vs['name'].index(synton1), graph.vs['name'].index(synton2))

def routine(graph, synton_of_targets, params, mode='A'):
    if mode == 'A':
        look_at = 0
    else:
        look_at = 1
    sorted_syntons = sorted(graph.vs['name'], key=lambda synton: synton[look_at])
    evalute_proximity(sorted_syntons, graph, params, look_at)
    cc = graph.components()
    [cc_with_target_target] = [cluster
                               for cluster in cc
                               if synton_of_targets in graph.vs[cluster]['name']]
    return cc_with_target_target

def compute_score(syntons, synton_of_targets, boolean_synton_of_targets):
    minGA = min([posA for posA, posB in syntons])
    maxGA = max([posA for posA, posB in syntons])
    minGB = min([posB for posA, posB in syntons])
    maxGB = max([posB for posA, posB in syntons])
    synt_coverage = (maxGA-minGA+1) + (maxGB-minGB+1)
    if not boolean_synton_of_targets:
        syntons.remove(synton_of_targets)
    geneA_in_synt = len(set([synton[0] for synton in syntons]))
    geneB_in_synt = len(set([synton[1] for synton in syntons]))
    avg_genes = (geneA_in_synt + geneB_in_synt)/2
    weighting = (geneA_in_synt + geneB_in_synt)/synt_coverage
    score = avg_genes * weighting
    return score

def find_common_connected_components(maxiG, gA, gB, targetA, targetB, synton_of_targets, boolean_synton_of_targets, params):
    ccA = routine(gA, synton_of_targets, params)
    ccB = routine(gB, synton_of_targets, params, mode='B')
    itrsect = list(set(ccA) & set(ccB))
    print(gB.vs[itrsect]['name'])
    if boolean_synton_of_targets:
        if len(itrsect) >= 2:
            print('are long enough !')
            print(ccA, ccB, itrsect)
            if ccA == ccB == itrsect:
                print('we have a ccc !!!')
                score = compute_score(gA.vs[itrsect]['name'], synton_of_targets, boolean_synton_of_targets)
                print(score)
                if not maxiG.vertex_attributes() or targetA not in maxiG.vs['name']:
                    maxiG.add_vertex(targetA)
                if not maxiG.vertex_attributes() or targetB not in maxiG.vs['name']:
                    maxiG.add_vertex(targetB)
                print(targetA, targetB)
                vertex_idx_targetA = maxiG.vs['name'].index(targetA)
                vertex_idx_targetB = maxiG.vs['name'].index(targetB)
                maxiG.add_edge(vertex_idx_targetA, vertex_idx_targetB)
                edge_idx_AB = maxiG.get_eid(vertex_idx_targetA, vertex_idx_targetB)
                maxiG.es[edge_idx_AB]['weight'] = score
            else:
                print('we need more steps')
                gA_memory = gA.copy()
                gA = ig.Graph()
                gA.add_vertices(gA_memory.vs[itrsect]['name'])
                gB = gA.copy()
                find_common_connected_components(maxiG, gA, gB, targetA, targetB, synton_of_targets, boolean_synton_of_targets, params)
        else:
            print('not long enough :\'{')

    else:
        if len(itrsect) > 2:
            print('are long enough !')
            print(ccA, ccB, itrsect)
            if ccA == ccB == itrsect:
                print('we have a ccc !!!')
                compute_score(itrsect, boolean_synton_of_targets)
                maxiG.add_vertex((targetA, targetB))
            else:
                print('we need more steps')
                gA_memory = gA.copy()
                gA = ig.Graph()
                gA.add_vertices(gA_memory.vs[itrsect]['name'])
                gB = gA.copy()
                find_common_connected_components(maxiG, gA, gB, targetA, targetB, synton_of_targets, boolean_synton_of_targets, params)
        else:
            print('not long enough :\'{')
    return maxiG

def build_maxi_graph(maxiG, targets_syntons, params):
    for (targetA, targetB) in targets_syntons:
        print('######################################')
        print(targetA, targetB)
        gA = ig.Graph()
        gA.add_vertices(targets_syntons[(targetA, targetB)]['syntons'])
        gB = gA.copy()
        synton_of_targets = targets_syntons[(targetA, targetB)]['AB_in_synt'][:-1]
        boolean_synton_of_targets = targets_syntons[(targetA, targetB)]['AB_in_synt'][-1]
        print('first search')
        maxiG = find_common_connected_components(maxiG, gA, gB, targetA, targetB, synton_of_targets, boolean_synton_of_targets, params)
        if maxiG.vertex_attributes():
            print(maxiG.vs['name'])
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
    with open('{}/test/TMP/{}/{}'.format(TMPDIRECTORY, 'ClusteringIntoFamilies', 'targets_list'), 'rb') as file:
        targets_list = pickle.load(file)

    targets_syntons = {}

    cds_info = set_userGC_similarityContext(targets_list, cds_info, params)
    for idx, targetA in enumerate(targets_list[:-1]):
        for targetB in targets_list[idx+1:]:
            if cds_info[targetA]['genome'] != cds_info[targetB]['genome']:
                families_intersect = intersect_families(cds_info[targetA],
                                                        cds_info[targetB])
                if len(families_intersect) > 1:
                    for afam in families_intersect:
                        A_in_synt = get_pos_synton(targetA, afam, cds_info)
                        B_in_synt = get_pos_synton(targetB, afam, cds_info)
                        syntons = cartesian_product(A_in_synt, B_in_synt)
                        #print(targetA, targetB, afam, syntons)
                        targets_syntons.setdefault((targetA, targetB),
                                                   {}).setdefault('syntons',
                                                                  []).extend(syntons)
                    add_target_target_synton(targetA, targetB, cds_info, targets_syntons)

    maxi_graph = ig.Graph()
    maxi_graph = build_maxi_graph(maxi_graph, targets_syntons, params)
    print(maxi_graph)
    print('\n')

    ### Edge-betweenness clustering
    # graph_edge_btwness = maxi_graph.community_edge_betweenness(directed=False)
    # print(graph_edge_btwness)
    # complete_dendogram = fix_dendrogram(maxi_graph, graph_edge_btwness)    
    # btwness_clusters = complete_dendogram.as_clustering()
    # print(btwness_clusters)
    
    ### Walktrap clustering
    print('\n*** Walktrap clustering ***')
    graph_walktrap = maxi_graph.community_walktrap(weights=maxi_graph.es['weight'])
    print(graph_walktrap)
    walktrap_clustering = graph_walktrap.as_clustering()
    print(walktrap_clustering)

    ### Louvain clustering
    print('\n*** Louvain clustering ***')
    graph_louvain = maxi_graph.community_multilevel()
    print(graph_louvain)


    #layout = maxi_graph.layout('circle')
    #layout = maxi_graph.layout('drl')
    #layout = maxi_graph.layout('large')
    #layout = maxi_graph.layout('random')
    #layout = maxi_graph.layout('tree')
    #layout = maxi_graph.layout('circular')

    # layout = maxi_graph.layout('fr') # pourquoi pas ...
    # #layout = maxi_graph.layout('kk') # pourquoi pas ...
    # ig.plot(maxi_graph, 'test_graph_fr.png', layout=layout, edge_width=maxi_graph.es['weight'], vertex_label=maxi_graph.vs['name'])


    #g1.add_vertices([target for duo_targets in targets_syntons for target in targets_syntons[duo_targets]])
    #print(g1)


    # target_famies = get_families_per_target(families, cds_info)
    # syntenies = intersect_target_famies(target_famies, cds_info, families)
    # synteny_filter(syntenies, params)
    # print('syntenies:\n{}'.format(syntenies))

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3], os.path.abspath('.'))
