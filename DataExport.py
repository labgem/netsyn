#!/usr/bin/env python3

##########
# Import #
##########
import argparse
#import json
#import pickle
import logging
import os
#import time
import igraph as ig
import common
#############
# Functions #
#############

def parse_tsv(fname, ac=None, mc=None):
    ''' create a list of dictionaries
        get from input file (fname) information by line
        every line is a dictionary
    '''
    logger = logging.getLogger('parse_tsv')
    first_line = True
    errors = False
    rows = []
    #replicons_index = {}
    #line_number = 0
    with open(fname, 'r') as file:
        for line in file:
            if first_line:
                headers = {}
                for index, header in enumerate(line.split('\t')):
                    header = header.replace('\r\n', '').replace('\n', '') # header.strip() ???
                    # p = re.compile(r'(?:{})'.format('|'.join(ac)))
                    # if not p.search(header):
                    #     logger.error('{}: Column name not valid.'.format(header))
                    #     errors = True
                    if header in headers.values():
                        logger.error('{}: Duplicated column.'.format(header))
                        errors = True
                    headers[index] = header
                # if mc:
                #     errors = check_headers(errors, mc, headers)
                first_line = False
            else:
                #line_number += 1
                row = {}
                #row['line_number'] = line_number
                for index, column in enumerate(line.split('\t')):
                    row[headers[index]] = column.replace('\r\n', '').replace('\n', '') # header.strip()
                rows.append(row)
    if errors:
        logger.error('Madatory columns')
        exit(1)
    return rows

def run(NODES, EDGES, TAXONOMY, CONTIGS, METADATA, RESULTSDIR):
    # Constants
    boxName = common.global_dict['boxName']['DataExport']
    # Outputs
    graphML = common.global_dict['files']['DataExport']['graphML']
    nodesRes = common.global_dict['files']['DataExport']['nodes']
    edgesRes = common.global_dict['files']['DataExport']['edges']
    #htmlOut = common.global_dict['files']['DataExport']['html']
    #settingsOut = common.global_dict['files']['DataExport']['settings']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    logger.info('{} running...'.format(boxName))
    # Process
    if not os.path.isdir(RESULTSDIR):
        os.mkdir(RESULTSDIR)

    list_of_nodes = common.read_pickle(NODES)
    list_of_edges = common.read_pickle(EDGES)
    taxonomicLineage = common.read_pickle(TAXONOMY)
    contigs = common.read_pickle(CONTIGS)
    metadata = parse_tsv(METADATA) # liste de dicos

    # COM: addition of metadata and taxonomic lineage information to the nodes
    g = ig.Graph()
    for idx, anode in enumerate(list_of_nodes):
        uniP = anode['UniProtAC']
        g.add_vertex(uniP)
        g.vs[idx]['inc_cds'] = anode['cds_ref']
        for adata in metadata:
            if adata['UniProtAC'] == uniP:
                anode.setdefault('MetaData', {})
                for key, value in adata.items():
                    if key != 'UniProtAC':
                        anode['MetaData'].update({key:value})
                        g.vs[idx][key] = value
                break
        contig_ref = anode['Contig']
        taxon_id = contigs[contig_ref]['taxon_ID']
        anode['Lineage'] = taxonomicLineage[taxon_id].copy()
        for key, value in anode['Lineage'].items():
            g.vs[idx][key] = value[1]
        for key, value in anode['Clustering'].items():
            g.vs[idx][key] = value

    for idx, aedge in enumerate(list_of_edges):
        source = aedge['source']
        target = aedge['target']
        g.add_edge(g.vs['name'].index(source), g.vs['name'].index(target))
        g.es[idx]['Families'] = aedge['families']
        g.es[idx]['weight'] = aedge['score']

    g.write_graphml(graphML)
    common.write_json(list_of_nodes, nodesRes)
    common.write_json(list_of_edges, edgesRes)

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(usage='''DataExport options...''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-n', '--Nodes', type=str,
                        required=True, help='Path of the nodes file obtained from the SyntenyFinder part')
    group1.add_argument('-e', '--Edges', type=str,
                        required=True, help='Path of the edges file obtained from the SyntenyFinder part')
    group1.add_argument('-taxo', '--Taxonomy', type=str,
                        required=True, help='Path of the taxonomyLineage file from the ClusteringIntoFamilies part')
    group1.add_argument('-c', '--ContigsInfo', type=str,
                        required=True, help='Path of the contigs file from the ClusteringIntoFamilies part')
    group1.add_argument('-md', '--MetaData', type=str,
                        required=True, help='Path of the metadata file provided by the user')
    group1.add_argument('-o', '--OutputName', type=str, required=True,
                        help='Output name files.')

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
    boxName = common.global_dict['boxName']['DataExport']
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('graphML', '{}_Results.graphML'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('nodes', '{}_Results_nodes.json'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('edges', '{}_Results_edges.json'.format(args.OutputName))
    #######
    # Run #
    #######
    run(args.Nodes, args.Edges, args.Taxonomy, args.ContigsInfo, args.MetaData, '.')
