#!/usr/bin/env python3

##########
# Import #
##########
import argparse
import json
import pickle
import logging
import os
import igraph as ig
import time
#############
# Functions #
#############

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='''Description of the DataExport usage''',
                                     epilog='''All's well that ends well.''',
                                     usage='''DataExport options...''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-n', '--Nodes', type=str,
                        required=True, help='Path of the nodes file obtained from the SyntenyFinder part')
    parser.add_argument('-e', '--Edges', type=str,
                        required=True, help='Path of the edges file obtained from the SyntenyFinder part')
    parser.add_argument('-taxo', '--Taxonomy', type=str,
                        required=True, help='Path of the taxonomyLineage file from the ClusteringIntoFamilies part')
    parser.add_argument('-c', '--ContigsInfo', type=str,
                        required=True, help='Path of the contigs file from the ClusteringIntoFamilies part')
    parser.add_argument('-md', '--MetaData', type=str,
                        required=True, help='Path of the metadata file provided by the user')
    parser.add_argument('-pn', '--ProjectName', type=str, required=True,
                         help='The project name.')
    return parser

def parse_tsv(fname, ac=None, mc=None):
    ''' create a list of dictionaries
        get from input file (fname) information by line
        every line is a dictionary
    '''
    logger = logging.getLogger('parse_tsv')
    first_line = True
    errors = False
    rows = []
    replicons_index = {}
    line_number = 0
    with open(fname,'r') as file:
        for line in file:
            if first_line:
                headers = {}
                for index, header in enumerate(line.split('\t')):
                    header = header.replace('\r\n','').replace('\n','') # header.strip() ???
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
                    row[headers[index]] = column.replace('\r\n','').replace('\n','') # header.strip()
                rows.append(row)
    if errors:
        logger.error('Madatory columns')
        sys.exit(1)
    return rows

def run(BOXNAME, TMPDIRECTORY, NODES, EDGES, TAXONOMY, CONTIGS, METADATA):
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    logger.info('{} running...'.format(BOXNAME))
    TMPDIRECTORYPROCESS = '{}/{}'.format(TMPDIRECTORY, BOXNAME)
    if not os.path.isdir(TMPDIRECTORYPROCESS):
        os.mkdir(TMPDIRECTORYPROCESS)

    with open(NODES, 'rb') as file:
        list_of_nodes = pickle.load(file)
    with open(EDGES, 'rb') as file:
        list_of_edges = pickle.load(file)
    with open(TAXONOMY, 'rb') as file:
        taxonomicLineage = pickle.load(file)
    with open(CONTIGS, 'rb') as file:
        contigs = pickle.load(file)

    metadata = parse_tsv(METADATA) # liste de dicos

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
        taxon_id = contigs[contig_ref]['taxon_id']
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

    g.write_graphml('{}/test.graphml'.format(TMPDIRECTORYPROCESS))

    with open('{}/nodes.json'.format(TMPDIRECTORYPROCESS), 'w') as file:
        json.dump(list_of_nodes, file, indent=4)
    with open('{}/edges.json'.format(TMPDIRECTORYPROCESS), 'w') as file:
        json.dump(list_of_edges, file, indent=4)

if __name__ == '__main__':
    parser = argumentsParser()
    args = parser.parse_args()
    BOXNAME = 'DataExport'
    TMPDIRECTORY = '{}/TMP'.format(args.ProjectName)
    #os.mkdir(TMPDIRECTORY)
    #MAXGCSIZE = 11
    run(BOXNAME, TMPDIRECTORY, args.Nodes, args.Edges, args.Taxonomy, args.ContigsInfo, args.MetaData)#MAXGCSIZE, args.WindowSize, args.SyntenyGap, args.SyntenyScoreCuttoff)
