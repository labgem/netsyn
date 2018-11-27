##########
# Import #
##########
import argparse
import re
import sys
import random
import string
import pickle
import json
import os
import math
import shutil
import xml.etree.ElementTree as ET
import subprocess
import urllib3
import logging
from Bio import SeqIO
#############
# Functions #
#############

#('{}/TMP/GetINSDCFiles/inputClusteringIntoFamiliesStep.tsv'.format(args.ProjectName), args, TMPDIRECTORY)

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description = '''My Description. And what a lovely description it is. ''',
                                     epilog = '''All's well that ends well.''',
                                     usage = '''ClusteringIntoFamilies options...''',
                                     formatter_class = argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', type = str,
                        required = True, help = 'Path of the input obtained from the GetINSDCFiles part')
    parser.add_argument('-id', '--Ident', type = float,
                        default = 0.3, help='Sequence identity.\nDefault value: 0.3.')
    parser.add_argument('-mc', '--MinCoverage', type = float,
                        default = 0.8, help='Minimal coverage allowed.\nDefault value: 0.8.')
    parser.add_argument('-pn', '--ProjectName', type = str, required = True,
                        help = 'The project name.')
    return parser

def skip_duplicates(iterable, key=lambda x: x):
    ''' removes duplicates from a list keeping the order of the elements
    Uses a generator
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

def check_headers(errors, mandatory_columns, headers):
    ''' check if all mandatory columns are provided
    '''
    logger = logging.getLogger('check_headers')
    for mandatory_column in mandatory_columns:
        if not mandatory_column in headers.values():
            logger.error('{}: Missing column!'.format(mandatory_column))
            errors = True
    return errors

def parse_tsv(fname, authorized_columns, mandatory_columns):
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
                    p = re.compile(r'(?:{})'.format('|'.join(authorized_columns)))
                    if not p.search(header):
                        logger.error('{}: Column name not valid.'.format(header))
                        errors = True
                    if header in headers.values():
                        logger.error('{}: Duplicated column.'.format(header))
                        errors = True
                    headers[index] = header
                errors = check_headers(errors, mandatory_columns, headers)
                first_line = False
            else:
                line_number += 1
                row = {}
                row['line_number'] = line_number
                for index, column in enumerate(line.split('\t')):
                    row[headers[index]] = column.replace('\r\n','').replace('\n','') # header.strip()
                rows.append(row)
    if errors:
        logger.error('Madatory columns')
        sys.exit(1)
    return rows

def read_rows(rows):
    ''' print rows information by row
    '''
    logger = logging.getLogger('read_rows')
    for row in rows:
        logger.debug('{}\n'.format(row))
    return 0

def create_d_input(d_rows):
    ''' formatting information by filename
    '''
    logger = logging.getLogger('{}.{}'.format(create_d_input.__module__, create_d_input.__name__))
    d_input = {}

    for arow in d_rows: # skips the first line, with headers
        filename = arow['nucleic_File_Name']
        contig_id = arow['nucleic_AC']
        d_input.setdefault(
            filename, {}).setdefault(contig_id, {}).setdefault('target_list', [])
        
        if arow['protein_AC'] not in d_input[filename][contig_id]['target_list']:
            d_input[filename][contig_id]['target_list'].append(arow['protein_AC'])

        d_input[filename]['protein_AC_field'] = arow['protein_AC_field']
        d_input[filename]['nucleic_File_Format'] = arow['nucleic_File_Format']

        if 'taxon_id' in arow.keys():
            if not d_input[filename][contig_id]:
                d_input[filename][contig_id] = arow['taxon_id']
            else:
                if d_input[filename][contig_id] != arow['taxon_id']:
                    logger.error('The nucleaotide accession {} already refers to a taxon ID {}'.format(contig_id, d_input[filename][contig_id]['taxon_id']))
                    exit(1)
        logger.debug('{}\t{}'.format(filename, d_input[filename]))
    return d_input

def check_and_get_input(input):
    '''
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    -*-*-*- faire le check en même temps que la création du dico -*-*-*-
    -*-*-*- revoir la liste des checks à faire -*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    '''
    logger = logging.getLogger('{}.{}'.format(check_and_get_input.__module__, check_and_get_input.__name__))

    ## from JL
    authorized_columns = ['protein_AC', 'protein_AC_field', 'nucleic_AC', 'nucleic_File_Format', 'nucleic_File_Name', 'taxon_id']
    mandatory_columns = ['protein_AC', 'protein_AC_field', 'nucleic_AC', 'nucleic_File_Format', 'nucleic_File_Name']
    d_rows = parse_tsv(input, authorized_columns, mandatory_columns)
    #read_rows(d_rows)
    ## End from JL
    d_input = create_d_input(d_rows)
    return d_input

def get_from_dbxref(aFeature, dbref):
    ''' retrieves the value from the dbxref list
    '''
    logger = logging.getLogger('{}.{}'.format(get_from_dbxref.__module__, get_from_dbxref.__name__))
    if dbref == 'taxon':
        pattern = 'taxon:'
    elif dbref == 'MaGe':
        pattern = 'MaGe:'
    elif dbref == 'UniProt':
        pattern = 'UniProt*'

    result = 'NA'
    if aFeature.qualifiers.get('db_xref'):
        for aRef in aFeature.qualifiers.get('db_xref'):
            if re.match(pattern, aRef):
                result = aRef.split(r':')[-1]
    return result

def get_uniq_value(aFeature, ref):
    ''' makes sure that when there is a value, it is not a sequence of values
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*- in that case a value is a sequence of values, it has to   -*-*-*-
    -*-*-*- send an error, or a warning (need to know which index to  -*-*-*-
    -*-*-*- select)                                                   -*-*-*-
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    '''
    logger = logging.getLogger('{}.{}'.format(get_uniq_value.__module__, get_uniq_value.__name__))
    result = aFeature.qualifiers.get(ref)
    if not result:
        try:
            result = get_from_dbxref(aFeature, ref)
            return result
        except:
            return 'NA'
    if len(result) == 1:
        return result[0]
    else:
        logger.warning('there are several values instead of a uniq value: {}'.format(result))
        # exit(1)

def get_required_value(func, aFeature, *args):
    ''' function called if the value is mandatory (protein ID (ident),
    sequence)
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*- an error has to be implemented so that any mandatory value -*-*-*-
    -*-*-*- couldn't be None or 'NA'                                   -*-*-*-
    -*-*-*- maybe do a try/except evaluation
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    '''
    logger = logging.getLogger('{}.{}'.format(get_required_value.__module__, get_required_value.__name__))
    if not func(aFeature, *args) or func(aFeature, *args) == 'NA':
        logger.error('A required value is not provided: {}'.format([arg for arg in args]))
        # exit(1)
        return 1
    else:
        return func(aFeature, *args)

def search_taxonID(aFeature, given_taxon):
    ''' makes sure that the taxon ID is provided by the user or through the
    INSDC file
    '''
    logger = logging.getLogger('{}.{}'.format(search_taxonID.__module__, search_taxonID.__name__))
    if given_taxon:
        #print('taxon_id provided by the user')
        taxon_id = given_taxon
    else:
        #print('taxon_id must be in the file')
        taxon_id = get_required_value(get_uniq_value, aFeature, 'taxon')
    return taxon_id

def get_contig_info(aFeature, contig_content, given_taxon):
    ''' then get the relied information to a contig as organism, strain, size
    and taxon ID
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*- an error has to be implemented so that any cds couldn't be -*-*-*-
    -*-*-*- relied to a taxon ID                                       -*-*-*-
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    '''
    logger = logging.getLogger('{}.{}'.format(get_contig_info.__module__, get_contig_info.__name__))
    taxon_id = search_taxonID(aFeature, given_taxon)
    if taxon_id != 'NA':
        contig_content.update({
                'organism': get_uniq_value(aFeature, 'organism'),
                'strain': get_uniq_value(aFeature, 'strain'),
                'taxon_id': taxon_id,
                'size': [aFeature.location.start.real+1,
                         aFeature.location.end.real
                         ],
                'cds_to_keep': [],
                'window': []
                })
    else:
        logger.error('We have a problem')
        exit(1)
    return contig_content

def is_pseudogene(aFeature):
    ''' tests if a CDS is a pseudogene
    '''
    return ('pseudo' in aFeature.qualifiers)
# if 'pseudo' in aFeature.qualifiers:
#     return True
# else:
#     return False

def get_pseudo_id(params):
    ''' creates an identifier for pseudogenes
    '''
    params['INC_PSEUDO_REF'] += 1
    INC_PSEUDO_REF = params['INC_PSEUDO_REF']
    return ':'.join(['PSEUDO', str(INC_PSEUDO_REF)])

def det_frame(direction, startCDS, endGenome):
    ''' determines the frame of the CDS
    '''
    if direction == '1':
        window = ((int(startCDS)-1)%3)+1
        strand = ''.join(['+', str(window)])
    elif direction == '-1':
        window = ((int(endGenome)-int(startCDS))%3)+1
        strand = ''.join(['-', str(window)])
    return strand

def get_pseudo_info(aFeature, cds_info, contig_content, params):
    ''' adds to the window's list information on the pseudogene
    '''
    logger = logging.getLogger('{}.{}'.format(get_pseudo_info.__module__, get_pseudo_info.__name__))
    params['INC_CDS_REF'] += 1
    INC_CDS_REF = params['INC_CDS_REF']

    pseudo_id = get_pseudo_id(params)
    start, stop = aFeature.location.start.real+1, aFeature.location.end.real
    cds_info[INC_CDS_REF] = {
        'protein_id': pseudo_id,
        'position': [start, stop],
        'frame': det_frame(str(aFeature.location.strand), start, contig_content.get('size')[1]),
        'product': get_uniq_value(aFeature, 'product'),
        'sequence': get_uniq_value(aFeature, 'translation'),
        'target': '',
        'context': None,
        'contig': INC_CONTIG_REF,
        'genome': params['INC_FILE'],
        'similarityFamily': None
        }
    contig_content['window'].append(INC_CDS_REF)
    return cds_info, contig_content, params

def get_prot_info(aFeature, cds_info, contig_content, proteinField, params):
    ''' adds to the window's list information on the protein
    '''
    logger = logging.getLogger('{}.{}'.format(get_prot_info.__module__, get_prot_info.__name__))
    INC_CONTIG_REF = params['INC_CONTIG_REF']
    params['INC_CDS_REF'] += 1
    INC_CDS_REF = params['INC_CDS_REF']

    ident = get_required_value(get_uniq_value, aFeature, proteinField)
    if ident == 1:
        ident = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(15))
        logger.info('Protein {} has no field protein_id. A random name has been created: {}'.format(INC_CDS_REF, ident))
    start, stop = aFeature.location.start.real+1, aFeature.location.end.real
    ec_nums = (aFeature.qualifiers.get('EC_number') if aFeature.qualifiers.get('EC_number') else 'NA')
    cds_info[INC_CDS_REF] = {
        'protein_id': ident,
        'position': [start, stop],
        'frame': det_frame(str(aFeature.location.strand), start, contig_content.get('size')[1]),
        'product': get_uniq_value(aFeature, 'product'),
        'sequence': get_required_value(get_uniq_value, aFeature, 'translation'),
        'ec_number': ec_nums,
        'uniprot': get_from_dbxref(aFeature, 'UniProt'),
        'target': [],
        'context': None,
        'contig': INC_CONTIG_REF,
        'genome': params['INC_FILE'],
        'similarityFamily': None
        }
    contig_content['window'].append(INC_CDS_REF)
    return cds_info, contig_content, params

def set_target_to_gc(ref_target, gcList, cds_info):
    ''' add the target reference to members of its genomic context
    '''
    for gc_member in gcList:
        if ref_target not in cds_info[gc_member]['target']:
            cds_info[gc_member]['target'].append(ref_target)
    return cds_info

def is_target(cds, target_list, cds_info, contig_content, targets_storage, params):
    ''' tests if the tested CDS is referenced as a target
    If yes, copy the window's information to cds_keeper :
    information on target and genomic context are saved

    *** tested value used to be the CDS in the middle of the window ***
    '''
    logger = logging.getLogger('{}.{}'.format(is_target.__module__, is_target.__name__))
    if cds_info[cds]['protein_id'] in target_list:
        targets_storage.append(cds)
        contig_content['cds_to_keep'].extend(contig_content['window'])

        cds_info[cds]['target'].append(cds)
        cds_info[cds]['context'] = contig_content['window'].copy()
        cds_info[cds]['similarityContext'] = None

        cds_info = set_target_to_gc(cds, cds_info[cds]['context'], cds_info)
        params['INC_TARGET_LOADED'] += 1
        logger.debug('cds ({}/{})\t- contig_content ({}) -\twindow: {}'.format(cds, cds_info[cds]['protein_id'], contig_content['contig'], contig_content['window']))
        logger.debug('cds ({}/{})\t- contig_content ({}) -\tcds_to_keep: {}'.format(cds, cds_info[cds]['protein_id'], contig_content['contig'], contig_content['cds_to_keep']))
    return cds_info, contig_content, targets_storage, params

def is_kept(cds, cds_info, contig_content):
    ''' verifies if the tested CDS is in cds_to_keep list
    If not, all information on it are removed

    *** tested value is the first value of the window ***
    '''
    logger = logging.getLogger('{}.{}'.format(is_kept.__module__, is_kept.__name__))
    if contig_content['cds_to_keep']:
        if cds not in [ref for ref in contig_content['cds_to_keep']]:
            del cds_info[cds]
    else:
        del cds_info[cds]
    del contig_content['window'][0]
    return cds_info, contig_content

def sliding_window(cds_info, contig_content, target_list, beginContig, targets_storage, params):
    ''' makes the window going ahead through the CDSs in a contig
    '''
    logger = logging.getLogger('{}.{}'.format(sliding_window.__module__, sliding_window.__name__))
    MAX_GC = params['MAX_GC']
    HALF_SIZE_GC = math.floor(MAX_GC/2)
    window = contig_content['window']
    window_size = len(window)
    if window_size == MAX_GC:
        beginContig = False
        cds_info, contig_content, targets_storage, params = is_target(window[HALF_SIZE_GC], target_list, cds_info, contig_content, targets_storage, params)
        cds_info, contig_content = is_kept(window[0], cds_info, contig_content)
    # COM: window's size didn't reach the MAX_GC yet
    elif window_size > HALF_SIZE_GC and beginContig:
        cds_info, contig_content, targets_storage, params = is_target(window[window_size-(HALF_SIZE_GC+1)], target_list, cds_info, contig_content, targets_storage, params)
    # COM: the window is overrunning the contig
    elif window_size > HALF_SIZE_GC and not beginContig:
        cds_info, contig_content, targets_storage, params = is_target(window[HALF_SIZE_GC], target_list, cds_info, contig_content, targets_storage, params)
        cds_info, contig_content = is_kept(window[0], cds_info, contig_content)
    # COM: there is no more target at the end of the contig
    elif window_size <= HALF_SIZE_GC and not beginContig:
        cds_info, contig_content = is_kept(window[0], cds_info, contig_content)
    return cds_info, contig_content, beginContig, targets_storage, params

def parse_insdc(afile, d_infile, cds_info, contig_info, targets_storage, params):
    ''' gets information for every contig mentionned in the input
    and gets information for every cds that is contained in a (MAX_GC sized-)
    window centered on a target
    '''
    logger = logging.getLogger('{}.{}'.format(parse_insdc.__module__, parse_insdc.__name__))
    STOP_INC_CONTIG = params['INC_CONTIG_REF'] + len(d_infile)-2

    typeParsing = d_infile['nucleic_File_Format']
    fieldProteinID = d_infile['protein_AC_field']

    with open(afile, 'r') as insdcFile:
        seqRecordParsed = SeqIO.parse(insdcFile, typeParsing)
        for seqRecord in seqRecordParsed:
            if params['INC_CONTIG_REF'] >= STOP_INC_CONTIG:
                logger.debug('To many contigs to parse ({}). Stop on contig {}'.format(params['INC_CONTIG_REF'], seqRecord.id))
                break
            contig_name = seqRecord.id#.split(r'.')[0]
            if contig_name in d_infile:
                params['INC_CONTIG_REF'] += 1
                INC_CONTIG_REF = params['INC_CONTIG_REF']
                params['INC_TARGET_LOADED'] = 0
                beginContig = True
                contig_info[INC_CONTIG_REF] = {'contig': contig_name}
                for aFeature in seqRecord.features:
                    if aFeature.type == 'source':
                        given_taxon_id = 'taxon_id' in d_infile[contig_name].keys()
                        contig_info[INC_CONTIG_REF] = get_contig_info(
                            aFeature,
                            contig_info[INC_CONTIG_REF],
                            given_taxon_id
                            )
                    elif aFeature.type == 'CDS':
                        newCdsAdded = False
                        if params['INC_TARGET_LOADED'] >= len(d_infile[contig_name]['target_list']):
                            logger.debug('No need to search more targets. {} already found on this contig {}'.format(params['INC_TARGET_LOADED'],seqRecord.id))
                            break
                        if is_pseudogene(aFeature):
                            if params['PSEUDOGENE']:
                                cds_info, contig_info[INC_CONTIG_REF], params = get_pseudo_info(aFeature, cds_info, contig_info[INC_CONTIG_REF], params)
                                newCdsAdded = True
                            #else: # ???
                                #break # ???
                        else:
                            cds_info, contig_info[INC_CONTIG_REF], params = get_prot_info(aFeature, cds_info, contig_info[INC_CONTIG_REF], fieldProteinID, params)
                            # print(params, [contig_info[ref]['window']
                            #                for ref in contig_info])
                            newCdsAdded = True
                        if newCdsAdded:
                            cds_info, contig_info[INC_CONTIG_REF], beginContig, targets_storage, params = sliding_window(
                                cds_info,
                                contig_info[INC_CONTIG_REF],
                                d_infile[contig_name]['target_list'],
                                beginContig,
                                targets_storage,
                                params
                                )
                beginContig = False
                for i in range(min(params['MAX_GC']-1, len(contig_info[INC_CONTIG_REF]['window']))):# don't work with '_' instead of 'i'
                    cds_info, contig_info[INC_CONTIG_REF], beginContig, targets_storage, params = sliding_window(
                        cds_info,
                        contig_info[INC_CONTIG_REF],
                        d_infile[contig_name]['target_list'],
                        beginContig,
                        targets_storage,
                        params
                        )
                #logger.debug('proteins referenced: {}'.format(cds_info.keys()))
                logger.debug('proteins in cds_to_keep on this contig: {}'.format(contig_info[INC_CONTIG_REF]['cds_to_keep']))
                contig_info[INC_CONTIG_REF]['cds_to_keep'] = list(skip_duplicates(contig_info[INC_CONTIG_REF]['cds_to_keep']))

    return cds_info, contig_info, targets_storage, params

def parse_INSDC_files(d_input, cds_info, contig_info, params):
    ''' every INSDC file in d_input will be parsed
        information are stored in dictionaries, one relative to contigs, the other relative to cds
        found targets are listed
    '''
    logger = logging.getLogger('{}.{}'.format(parse_INSDC_files.__module__, parse_INSDC_files.__name__))
    targets_storage = []
    nbr_of_files = len(d_input)
    for afile in d_input:
        params['INC_FILE'] += 1
        logger.info('Parsing of INSDC file ({}/{}): {}'.format(params['INC_FILE'], nbr_of_files, afile))
        cds_info, contig_info, targets_storage, params = parse_insdc(afile, d_input[afile], cds_info, contig_info, targets_storage, params)
        #print(params['INC_CONTIG_REF'], len(d_input[afile])-2)
    return cds_info, contig_info, targets_storage, params

def concat_by_dot(alist):
    ''' does the concatenation by a dot
    '''
    return '.'.join(alist)

def write_multiFasta(cds_info, output):
    ''' writes a multiFasta file from the dictionary obtained by the INSDC file
    parser
    '''
    logger = logging.getLogger('{}.{}'.format(write_multiFasta.__module__, write_multiFasta.__name__))
    #fileToWrite = concat_by_dot(prefix, suffix)
    with open(output, 'w') as fastaFile:
        for cds in cds_info:
            count = 0
            fastaFile.write('>{}\n'.format(cds))
            for ac_amine in cds_info[cds]['sequence']:
                count += 1
                fastaFile.write(ac_amine)
                if count >= 80:
                    fastaFile.write('\n')
                    count = 0
            fastaFile.write('\n')
    # if os.path.getsize(output) == 0:
    #     logger.debug('MultiFasta File is empty.')
    #     exit(1)
    return 0

def write_pickle(dictionary, output):
    ''' writes all dictionaries of INSDC file parsed in a pickle file format
    '''
    logger = logging.getLogger('{}.{}'.format(write_pickle.__module__, write_pickle.__name__))
    with open(output, 'wb') as pickleFile:
        pickle.dump(dictionary, pickleFile)
    return 0

def write_json(dictionary, output):
    ''' writes all dictionaries of INSDC file parsed in a json file format
    '''
    logger = logging.getLogger('{}.{}'.format(write_json.__module__, write_json.__name__))
    with open(output, 'w') as jsonFile:
        json.dump(dictionary, jsonFile, indent=4)
    return 0

# def write_list(list, output):
#     logger = logging.getLogger('{}.{}'.format(write_list.__module__, write_list.__name__))
#     with open(output, 'w') as opt:
#         for atarget in list:
#             opt.write('{}\n'.format(atarget))
#     return 0

def get_lineage(xml):
    ''' extracts the taxonomic lineage from the provided xml file format
    '''
    logger = logging.getLogger('{}.{}'.format(get_lineage.__module__, get_lineage.__name__))
    lineage_full = {}
    root = ET.fromstring(xml)
    scientificName = root.find('taxon').get('scientificName')
    rank = False
    if root.find('taxon').get('rank'):
        rank = root.find('taxon').get('rank')
    taxId = root.find('taxon').get('taxId')
    lineage_full[scientificName] = {'rank' : rank, 'taxId' : taxId}

    lineage = root.find('taxon').find('lineage')
    taxons = lineage.findall('taxon')
    for taxon in taxons:
        scientificName = taxon.get('scientificName')
        taxId = taxon.get('taxId')
        rank = False
        if taxon.get('rank'):
            rank = taxon.get('rank')
        lineage_full[scientificName] = {'rank' : rank, 'taxId' : taxId}
    return lineage_full

def get_taxo_from_web(taxonID):
    ''' does web request to get the taxonomic lineage for a taxon ID
    '''
    logger = logging.getLogger('{}.{}'.format(get_taxo_from_web.__module__, get_taxo_from_web.__name__))
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    http = urllib3.PoolManager()
    xml = http.request('GET', 'https://www.ebi.ac.uk/ena/data/view/Taxon:' + str(taxonID) + '&display=xml')
    with open('xml_file', 'w') as file:
        file.write(xml.data.decode('utf-8'))
    taxonLineage = get_lineage(xml.data.decode('utf-8'))
    return taxonLineage

def get_desired_lineage(lineage_full):
    ''' takes the ranks of interest through the taxonomic lineage
    '''
    desired_ranks = {
        'superkingdom' : 1,
        'phylum' : 5,
        'class' : 8,
        'order' : 13,
        'family' : 17,
        'genus' : 21,
        'species' : 26
    }
    desired_lineage = []

    for scientificName in lineage_full:
        if lineage_full[scientificName]['rank'] in desired_ranks:
            rank = lineage_full[scientificName]['rank']
            level = desired_ranks[rank]
            taxId = lineage_full[scientificName]['taxId']
            desired_lineage.append([level, rank, scientificName, taxId])
            del desired_ranks[rank]
    if desired_ranks != {}:
        for rank in desired_ranks:
            level = desired_ranks[rank]
            scientificName = 'NA'
            taxId = 'NA'
            desired_lineage.append([level, rank, scientificName, taxId])
    return desired_lineage

def store_into_dict(ataxon, desired_taxo):
    ''' puts information on taxonomic lineage into a dictionary
    '''
    d_newLineage = {}
    d_newLineage[ataxon] = {}
    for alevel in desired_taxo:
        d_newLineage[ataxon][alevel[1]] = [level_info
                                           for level_info in alevel
                                           if alevel.index(level_info) != 1]
    return d_newLineage

def get_taxonLineage(taxonIDs):
    ''' concatenates the various taxonomic lineage dictionaries in a
    super-dictionary called all_lineages
    '''
    all_lineages = {}
    print(taxonIDs)
    for ataxon in taxonIDs:
        res_taxo = get_taxo_from_web(ataxon)
        desired_taxo = get_desired_lineage(res_taxo)
        new_lineage = store_into_dict(ataxon, desired_taxo)
        all_lineages.update(new_lineage) # no risk of overwriting
    return all_lineages

def mmseqs_createdb(TMPDIRECTORYPROCESS, prefix):
    ''' creates a database using the mmseqs software
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_createdb.__module__, mmseqs_createdb.__name__))
    suffix = 'faa'
    multiFasta = '{}/{}'.format(TMPDIRECTORYPROCESS, concat_by_dot([prefix, suffix]))
    with open('{}/{}'.format(TMPDIRECTORYPROCESS, 'mmseqs_createdb.log'), 'w') as file:
        db_creation = subprocess.run(['mmseqs', 'createdb', multiFasta, '{}/{}.{}'.format(TMPDIRECTORYPROCESS, prefix, 'DB')], stdout=file, stderr=file, check=True)
        logger.info('createdb - exit code: {}'.format(db_creation.returncode))
    return 0

def mmseqs_clustering(TMPDIRECTORYPROCESS, prefix, cov, ident, cov_mode):
    ''' does the clustering using the mmseqs software
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_clustering.__module__, mmseqs_clustering.__name__))
    suffixDB = 'DB'
    suffixCluster = 'cluster'
    if os.path.isdir('{}/{}'.format(TMPDIRECTORYPROCESS, 'MMseqsTMP')):
        shutil.rmtree('{}/{}'.format(TMPDIRECTORYPROCESS, 'MMseqsTMP'))
    os.mkdir('{}/{}'.format(TMPDIRECTORYPROCESS, 'MMseqsTMP'))
    dataBase = '{}/{}'.format(TMPDIRECTORYPROCESS, concat_by_dot([prefix, suffixDB]))
    outputCluster = '{}/{}'.format(TMPDIRECTORYPROCESS, concat_by_dot([prefix, suffixCluster]))
    with open('{}/{}'.format(TMPDIRECTORYPROCESS, 'mmseqs_clustering.log'), 'w') as file:
        clust_creation = subprocess.run(['mmseqs', 'cluster', dataBase,
                                         outputCluster, '{}/{}'.format(TMPDIRECTORYPROCESS, 'MMseqsTMP'),
                                         '--min-seq-id', str(ident),
                                         '--cov-mode', str(cov_mode),
                                         '-c', str(cov),
                                         '--cluster-mode', str(0)#str(clust_mode)#,
                                         #'--cascaded', str(cascaded)
                                         ], stdout=file, stderr=file, check=True)
        logger.info('clustering - exit code: {}'.format(clust_creation.returncode))
    return 0

def mmseqs_createTSV(TMPDIRECTORYPROCESS, prefix):
    ''' executes the mmseqs command line 'mmseqs createtsv'
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_createTSV.__module__, mmseqs_createTSV.__name__))
    suffixDB = 'DB'
    suffixCluster = 'cluster'
    suffixTSV = 'tsv'
    inputDB = '{}/{}'.format(TMPDIRECTORYPROCESS, concat_by_dot([prefix, suffixDB]))
    inputCluster = '{}/{}'.format(TMPDIRECTORYPROCESS, concat_by_dot([prefix, suffixCluster]))
    outputTSV = '{}/{}'.format(TMPDIRECTORYPROCESS, concat_by_dot([prefix, suffixTSV]))
    with open('{}/{}'.format(TMPDIRECTORYPROCESS, 'mmseqs_createtsv.log'), 'w') as file:
        tsv_creation = subprocess.run(['mmseqs', 'createtsv', inputDB,
                                       inputDB, inputCluster, outputTSV
                                       ], stdout=file, stderr=file, check=True)
        logger.info('createTSV - exit code: {}'.format(tsv_creation.returncode))
    return 0

def mmseqs_runner(params, TMPDIRECTORYPROCESS):
    ''' runs the mmseqs2 software on the multiFasta file 'ALL.faa'
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_runner.__module__, mmseqs_runner.__name__))
    logger.info('MMseqs2 running ...')
    mmseqs_createdb(TMPDIRECTORYPROCESS, params['prefix'])
    mmseqs_clustering(TMPDIRECTORYPROCESS, params['prefix'], params['coverage'],
                      params['min_id'], params['cov_mode'])
    mmseqs_createTSV(TMPDIRECTORYPROCESS, params['prefix'])
    shutil.rmtree('{}/{}'.format(TMPDIRECTORYPROCESS, 'MMseqsTMP/'))
    logger.info('End of MMseqs2 running !')
    return 0

def regroup_families(tsv_file, cds_info):
    ''' creates a dictionary to store families obtained by MMseqs2
    '''
    logger = logging.getLogger('{}.{}'.format(regroup_families.__module__, regroup_families.__name__))
    INC_FAMILY = 1
    centroid = None
    with open(tsv_file, 'r') as file:
        lines = file.readlines()
    for aline in lines:
        aline = aline.strip().split('\t')
        if not centroid:
            centroid = aline[0]
        elif centroid != aline[0]:
            centroid = aline[0]
            INC_FAMILY += 1
        cds = int(aline[1])
        cds_info[cds]['similarityFamily'] = INC_FAMILY 
    return cds_info

def run(BOXNAME, TMPDIRECTORY, INPUT_II, MAXGCSIZE, IDENT, COVERAGE):
    ''' main script to run the second box of NetSyn2
    '''
    #BOXNAME = 'ClusteringIntoFamilies'
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    logger.info('{} running...'.format(BOXNAME))
    TMPDIRECTORYPROCESS = '{}/{}'.format(TMPDIRECTORY, BOXNAME)
    if not os.path.isdir(TMPDIRECTORYPROCESS):
        os.mkdir(TMPDIRECTORYPROCESS)

    params = {
        'PSEUDOGENE': False, # Tells if pseudogenes are included in the analysis
        'MAX_GC': MAXGCSIZE, # size of the window
        'INC_PSEUDO_REF': 0, # counter of pseusogenes
        'INC_CDS_REF': 0,
        'INC_CONTIG_REF': 0,
        'INC_FILE': 0,
        'min_id': IDENT,
        'cov_mode': 1,
        'coverage': COVERAGE,
        }
    params['prefix'] = "MMseqs2_run"

    try:
        os.remove('{}/{}.{}'.format(TMPDIRECTORYPROCESS, params['prefix'], 'cluster'))
    except:
        pass
    try:
        os.remove('{}/{}.{}'.format(TMPDIRECTORYPROCESS, params['prefix'], 'cluster.index'))
    except:
        pass
    try:
        os.remove('{}/{}.{}'.format(TMPDIRECTORYPROCESS, params['prefix'], 'tsv'))
    except:
        pass
    try:
        os.remove('{}/{}'.format(TMPDIRECTORYPROCESS, 'mmseqs_createtsv.log'))
    except:
        pass

    cds_info = {}
    contig_info = {}
    d_input = check_and_get_input(INPUT_II)
    #tester la fonction map() de python pour appliquer une fonction sur une
    #liste
    #usage : map(myFun, myList)
    logger.info('INSDC files parsing ...')
    cds_info, contig_info, targets_storage, params = parse_INSDC_files(d_input, cds_info, contig_info, params)
    logger.info('End of INSDC files parsing !')
    # for akey in cds_info:
    #     print(akey, cds_info[akey], sep='\t')
    #     print('***')
    # print('\n~~~~~~~~~~~~~~\n')
    # for akey in contig_info:
    #     print(akey, contig_info[akey], sep='\t')
    #     print('***')

    logger.info('will write the multifasta file')
    write_multiFasta(cds_info, '{}/{}'.format(TMPDIRECTORYPROCESS, concat_by_dot([params["prefix"], 'faa'])))
    write_pickle(contig_info, '{}/{}'.format(TMPDIRECTORYPROCESS, 'contigs.pickle'))
    write_pickle(targets_storage, '{}/{}'.format(TMPDIRECTORYPROCESS, 'targets_list.pickle'))
    write_json(contig_info, '{}/{}'.format(TMPDIRECTORYPROCESS, 'contigs.json'))
    write_json(targets_storage, '{}/{}'.format(TMPDIRECTORYPROCESS, 'targets_list.json'))

    logger.debug('list of targets: {}'.format(targets_storage))
    logger.info('Written files:\n{}\n{}\n{}'.format(concat_by_dot([params["prefix"], 'faa']), 'genomicContexts.pickle', 'contigs.pickle'))

    #print(contig_info)
    taxonIDs = list(set([contig_info[contig]['taxon_id'] for contig in contig_info]))
    taxonomicLineage = get_taxonLineage(taxonIDs)
    write_pickle(taxonomicLineage, '{}/{}'.format(TMPDIRECTORYPROCESS, 'taxonomyLineage.pickle'))
    write_json(taxonomicLineage, '{}/{}'.format(TMPDIRECTORYPROCESS, 'taxonomyLineage.json'))
    logger.info('Written file:\n{}'.format('taxonomyLineage.pickle'))

    mmseqs_runner(params, TMPDIRECTORYPROCESS)
    
    cds_info = regroup_families('{}/{}'.format(TMPDIRECTORYPROCESS, concat_by_dot([params["prefix"], 'tsv'])), cds_info)
    write_pickle(cds_info, '{}/{}'.format(TMPDIRECTORYPROCESS, 'genomicContexts.pickle'))
    write_json(cds_info, '{}/{}'.format(TMPDIRECTORYPROCESS, 'genomicContexts.json'))
    # real_families = {key: value for (key, value) in families.items() if len(value) > 1}
    # singletons = {key: value for (key, value) in families.items() if len(value) == 1}
    # write_pickle(real_families, '{}/{}'.format(TMPDIRECTORYPROCESS, 'proteinFamilies.pickle'))
    # write_pickle(singletons, '{}/{}'.format(TMPDIRECTORYPROCESS, 'proteinSingletons.pickle'))
    # logger.info('Written files:\n{}\n{}'.format('proteinFamilies.pickle', 'proteinSingletons.pickle'))
    
    with open('{}/analysed_targets'.format(TMPDIRECTORY), 'w') as file:
        for inc in cds_info.keys():
            file.write('{}\t{}\t{}\n'.format(inc, cds_info[inc]['protein_id'], cds_info[inc]['uniprot']))

    print('Nbr of targets: {} - {}'.format(len(targets_storage), params['INC_TARGET_LOADED']))
    with open('analysed_targets', 'w') as file:
        file.write('List of targets:\n')
        for tar in targets_storage:
            file.write('{}\t{}\t{}\n'.format(tar, cds_info[tar]['protein_id'], cds_info[tar]['uniprot']))
    print('Nbr of cds: {}'.format(len(cds_info)))
    print('Nbr of parsed files: {}'.format(params['INC_FILE']))


    logger.info('End of ClusteringIntoFamilies')

if __name__ == '__main__':
    parser = argumentsParser()
    args = parser.parse_args()
    #def run(BOXNAME, TMPDIRECTORY, INPUT_II, MAXGCSIZE, IDENT, COVERAGE):
    print(args)
    BOXNAME = 'ClusteringIntoFamilies'
    TMPDIRECTORY = '{}/TMP'.format(args.ProjectName)
    MAXGCSIZE = 11
    run(BOXNAME, TMPDIRECTORY, args.input, MAXGCSIZE, args.Ident, args.MinCoverage)
