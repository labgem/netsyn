##########
# Import #
##########
import common
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
    parser = argparse.ArgumentParser(description='''My Description. And what a lovely description it is. ''',
                                     epilog='''All's well that ends well.''',
                                     usage='''ClusteringIntoFamilies options...''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', type=str,
                        required=True, help='Path of the input obtained from the GetINSDCFiles part')
    parser.add_argument('-id', '--Ident', type=float,
                        default=0.3, help='Sequence identity.\nDefault value: 0.3.')
    parser.add_argument('-mc', '--MinCoverage', type=float,
                        default=0.8, help='Minimal coverage allowed.\nDefault value: 0.8.')
    parser.add_argument('-pn', '--ProjectName', type=str, required=True,
                        help='The project name.')
    return parser

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
        every line is a dictionary stored in a list
    '''
    logger = logging.getLogger('parse_tsv')
    first_line = True
    errors = False
    rows = []
    line_number = 0
    with open(fname, 'r') as file:
        for line in file:
            if first_line:
                headers = {}
                for index, header in enumerate(line.split('\t')):
                    header = header.replace('\r\n', '').replace('\n', '') # header.strip() ???
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
                    row[headers[index]] = column.replace('\r\n', '').replace('\n', '') # header.strip()
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
        input: list of dictionaries
        output: dictionary 
    '''
    logger = logging.getLogger('{}.{}'.format(create_d_input.__module__, create_d_input.__name__))
    d_input = {}
    errors = False

    for arow in d_rows: # skips the first line, with headers
        filename = arow['nucleic_File_Name']
        contig_id = arow['nucleic_AC']
        d_input.setdefault(
            filename, {}).setdefault(contig_id, {}).setdefault('target_list', [])

        if arow['protein_AC'] not in d_input[filename][contig_id]['target_list']:
            d_input[filename][contig_id]['target_list'].append(arow['protein_AC'])

        d_input[filename]['protein_AC_field'] = arow['protein_AC_field']
        d_input[filename]['nucleic_File_Format'] = arow['nucleic_File_Format']

        if 'taxon_ID' in arow.keys() and arow['taxon_ID'] != common.global_dict['defaultValue']:
            if 'taxon_ID' not in d_input[filename][contig_id]:
                d_input[filename][contig_id]['taxon_ID'] = arow['taxon_ID']
            else:
                if d_input[filename][contig_id]['taxon_ID'] != arow['taxon_ID']:
                    logger.error('Taxon ID provided in line {} is different than a previously provided one for the same INSDC file {}'.format(arow['line_number'], filename))
                    errors = True
    if errors:
        logger.info('Taxon ID inconsistency. Please review your input file')
        exit(1)
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
    authorized_columns = common.global_dict['inputIIheaders'] + ['taxon_ID']
    mandatory_columns = ['protein_AC', 'protein_AC_field', 'nucleic_AC', 'nucleic_File_Format', 'nucleic_File_Name']
    d_rows = parse_tsv(input, authorized_columns, mandatory_columns)
    #read_rows(d_rows)
    ## End from JL
    d_input = create_d_input(d_rows)
    return d_input

def get_from_dbxref(aFeature, dbref):
    ''' retrieve the value from the dbxref list
    '''
    logger = logging.getLogger('{}.{}'.format(get_from_dbxref.__module__, get_from_dbxref.__name__))
    if dbref == 'taxon':
        pattern = 'taxon:'
    elif dbref == 'MaGe':
        pattern = 'MaGe:'
    elif dbref == 'UniProt':
        pattern = 'UniProt*'

    result = common.global_dict['defaultValue']
    if aFeature.qualifiers.get('db_xref'):
        for aRef in aFeature.qualifiers.get('db_xref'):
            if re.match(pattern, aRef):
                result = aRef.split(r':')[-1]
    return result

def get_uniq_value(aFeature, ref):
    ''' make sure that when there is a value, it is not a sequence of values
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
            return common.global_dict['defaultValue']
    if len(result) == 1:
        return result[0]
    else:
        logger.warning('there are several values instead of a uniq value for {} field: {}'.format(ref, result))
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
    if not func(aFeature, *args) or func(aFeature, *args) == common.global_dict['defaultValue']:
        logger.error('A required value is not provided: {}'.format([arg for arg in args]))
        # exit(1)
        return 1
    else:
        return func(aFeature, *args)

def get_taxonID(aFeature):
    ''' get taxon ID provided in a INSDC file
    '''
    logger = logging.getLogger('{}.{}'.format(get_taxonID.__module__, get_taxonID.__name__))
    taxon_ID = get_uniq_value(aFeature, 'taxon')
    return taxon_ID

def get_contig_info(aFeature, contig_content, taxon_ID):
    ''' get the relied information to a contig as organism, strain, size
    and taxon ID
    input: source feature from an INSDC file
    output: updated contig dictionary
    '''
    logger = logging.getLogger('{}.{}'.format(get_contig_info.__module__, get_contig_info.__name__))
    if taxon_ID == common.global_dict['defaultValue']:
        taxon_ID = get_taxonID(aFeature)
    contig_content.update({
            'organism': get_uniq_value(aFeature, 'organism'),
            'strain': get_uniq_value(aFeature, 'strain'),
            'taxon_ID': taxon_ID,
            'size': [aFeature.location.start.real+1,
                     aFeature.location.end.real
                     ],
            'cds_to_keep': [],
            'window': []
            })
    return contig_content

def is_pseudogene(aFeature):
    ''' test if a CDS is a pseudogene
    '''
    return ('pseudo' in aFeature.qualifiers)
# if 'pseudo' in aFeature.qualifiers:
#     return True
# else:
#     return False

def get_pseudo_id(params):
    ''' create an identifier for pseudogenes
    '''
    params['INC_PSEUDO_REF'] += 1
    INC_PSEUDO_REF = params['INC_PSEUDO_REF']
    return ':'.join(['PSEUDO', str(INC_PSEUDO_REF)])

def det_frame(direction, startCDS, endGenome):
    ''' determine the frame of the CDS
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*- needs to be rethought because of cases where there is a   -*-*-*-
    -*-*-*- join(##..##, ##..###)                                     -*-*-*-
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    '''
    if direction == '1':
        window = ((int(startCDS)-1)%3)+1
        strand = ''.join(['+', str(window)])
    elif direction == '-1':
        window = ((int(endGenome)-int(startCDS))%3)+1
        strand = ''.join(['-', str(window)])
    return strand

def get_pseudo_info(aFeature, cds_info, contig_content, params):
    ''' add information of the pseudogene to cds_info dictionary
    input: CDS feature having a '/pseudo' qualifier
    output: updated cds_info dictionary, window list of the contig_content
    modified to contain afeature
    '''
    logger = logging.getLogger('{}.{}'.format(get_pseudo_info.__module__, get_pseudo_info.__name__))
    params['INC_CDS_REF'] += 1
    INC_CDS_REF = params['INC_CDS_REF']
    INC_CONTIG_REF = params['INC_CONTIG_REF']

    pseudo_id = get_pseudo_id(params)
    start, stop = aFeature.location.start.real+1, aFeature.location.end.real
    cds_info[INC_CDS_REF] = {
        'protein_id': pseudo_id,
        'position': [start, stop],
        'frame': det_frame(str(aFeature.location.strand), start, contig_content.get('size')[1]),
        'product': aFeature.qualifiers.get('product') if aFeature.qualifiers.get('product') else common.global_dict['defaultValue'],
        'sequence': get_uniq_value(aFeature, 'translation'),
        'target': [],
        'contig': INC_CONTIG_REF,
        'genome': params['INC_FILE'],
        }
    contig_content['window'].append(INC_CDS_REF)
    return cds_info, contig_content, params

def get_prot_info(aFeature, cds_info, contig_content, proteinField, params):
    ''' add information of the protein to cds_info dictionary
    input: CDS feature of an INSDC file
    output: updated cds_info dictionary, window list of the contig_content
    modified to contain afeature
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
    ec_nums = (aFeature.qualifiers.get('EC_number') if aFeature.qualifiers.get('EC_number') else common.global_dict['defaultValue'])
    cds_info[INC_CDS_REF] = {
        'protein_id': ident,
        'position': [start, stop],
        'frame': det_frame(str(aFeature.location.strand), start, contig_content.get('size')[1]),
        'product': aFeature.qualifiers.get('product') if aFeature.qualifiers.get('product') else common.global_dict['defaultValue'],
        'sequence': get_required_value(get_uniq_value, aFeature, 'translation'),
        'ec_number': ec_nums,
        'uniprot': get_from_dbxref(aFeature, 'UniProt'),
        'target': [],
        'contig': INC_CONTIG_REF,
        'genome': params['INC_FILE'],
        }
    contig_content['window'].append(INC_CDS_REF)
    return cds_info, contig_content, params

def set_target_to_gc(ref_target, gcList, cds_info):
    ''' add the target reference to members of its genomic context
    input: target reference, list of target's genomic context
    output: updated cds_info dictionary on genomic context cds key
    '''
    for gc_member in gcList:
        #if ref_target not in cds_info[gc_member]['target']:
        cds_info[gc_member]['target'].append(ref_target)
    return cds_info

def found_target_procedure(target, cds_info, contig_content):
    ''' add the genomic context to the 'cds_to_keep' list, copy the 'window'
    to the 'context'
    input: target reference
    output: updated cds_info on target's context and updated contig_content on
    cds_to_keep
    '''
    logger = logging.getLogger('{}.{}'.format(found_target_procedure.__module__, found_target_procedure.__name__))
    contig_content['cds_to_keep'].extend(contig_content['window'])

    cds_info[target]['context'] = contig_content['window'].copy()

    cds_info = set_target_to_gc(target, cds_info[target]['context'], cds_info)
    logger.debug('target ({}/{})\t- contig_content ({}) -\twindow: {}'.format(target, cds_info[target]['protein_id'], contig_content['contig'], contig_content['window']))
    logger.debug('target ({}/{})\t- contig_content ({}) -\tcds_to_keep: {}'.format(target, cds_info[target]['protein_id'], contig_content['contig'], contig_content['cds_to_keep']))
    return cds_info, contig_content

def remove_useless_cds(cds, cds_info, contig_content):
    ''' remove the cds if is not in cds_to_keep
    input: cds to test
    output: potential modified cds_info dictionary where cds entry is removed
    '''
    logger = logging.getLogger('{}.{}'.format(remove_useless_cds.__module__, remove_useless_cds.__name__))
    if cds not in [ref for ref in contig_content['cds_to_keep']]:
        del cds_info[cds]
    return cds_info

def parse_insdc(afile, d_infile, cds_info, contig_info, params):
    ''' get information for every contig mentionned in the input
    and get information for every cds that is contained in a (MAX_GC sized-)
    window centered on a target
    input: INSDC file, file input information in d_infile dictionary
    output: cds_info dictionary containing the cds belonging to the genomic context
    of the targets and targets itself, and
    contig_info dictionary containing the explored contig information
    (see get_contig_info function), and
    params is modified where values are incremented (INC_CONTIG_REF, INC_CDS_REF)
    '''
    logger = logging.getLogger('{}.{}'.format(parse_insdc.__module__, parse_insdc.__name__))
    
    NOT_CONTIG = ['protein_AC_field', 'nucleic_File_Format'] #### ~~~~~~~~ Probably should be changed ~~~~~~~~~~ #####
    CONTIG_LIST = [contig for contig in d_infile if contig not in NOT_CONTIG]
    logger.debug('Contig list: {}'.format(CONTIG_LIST))

    MAX_GC = params['MAX_GC']
    HALF_SIZE_GC = math.floor(MAX_GC/2)
    typeParsing = d_infile['nucleic_File_Format']
    fieldProteinID = d_infile['protein_AC_field']

    with open(afile, 'r') as insdcFile:
        seqRecordParsed = SeqIO.parse(insdcFile, typeParsing)
        CONTIG = next(seqRecordParsed)
        while CONTIG_LIST and CONTIG:
            # COM: getting the contig ID depending on its format
            if CONTIG.id in d_infile:
                contig_name = CONTIG.id
            elif CONTIG.id.split(r'.')[0] in d_infile:
                contig_name = CONTIG.id.split(r'.')[0]
            else:
                contig_name = ''

            if contig_name:
                # COM: initialization or update off incremented values
                CONTIG_LIST.remove(contig_name)
                params['INC_CONTIG_REF'] += 1
                INC_CONTIG_REF = params['INC_CONTIG_REF']
                params['INC_TARGET_LOADED'] = 0
                contig_info[INC_CONTIG_REF] = {'contig': contig_name}
                TARGET_LIST = d_infile[contig_name]['target_list']
                logger.debug('Target list for {} contig: {}'.format(contig_name, TARGET_LIST))

                # COM: beginning of the parsing
                for aFeature in CONTIG.features:
                    if aFeature.type == 'source':
                        contig_info[INC_CONTIG_REF] = get_contig_info(
                            aFeature,
                            contig_info[INC_CONTIG_REF],
                            d_infile[contig_name]['taxon_ID'] if 'taxon_ID' in d_infile[contig_name].keys() else common.global_dict['defaultValue']
                            )
                    elif aFeature.type == 'CDS':
                        newCdsAdded = False
                        if is_pseudogene(aFeature):
                            if params['PSEUDOGENE']:
                                cds_info, contig_info[INC_CONTIG_REF], params = get_pseudo_info(aFeature, cds_info, contig_info[INC_CONTIG_REF], params)
                                newCdsAdded = True
                        else:
                            cds_info, contig_info[INC_CONTIG_REF], params = get_prot_info(aFeature, cds_info, contig_info[INC_CONTIG_REF], fieldProteinID, params)
                            newCdsAdded = True

                        if newCdsAdded:
                            window = contig_info[INC_CONTIG_REF]['window']
                            window_length = len(window)
                            contig_content = contig_info[INC_CONTIG_REF]
                            if window_length == MAX_GC:
                                presumed_target = window[HALF_SIZE_GC]
                                presumed_kicked_out_cds = window[0]
                                if cds_info[presumed_target]['protein_id'] in TARGET_LIST:
                                    params['INC_TARGET_LOADED'] += 1
                                    TARGET_LIST.remove(cds_info[presumed_target]['protein_id'])
                                    cds_info, contig_content = found_target_procedure(presumed_target, cds_info, contig_content)
                                if contig_content['cds_to_keep']:
                                    cds_info = remove_useless_cds(presumed_kicked_out_cds, cds_info, contig_content)
                                else:
                                    del cds_info[presumed_kicked_out_cds]
                                del window[0]
                            elif window_length > HALF_SIZE_GC:
                                presumed_target = window[window_length-(HALF_SIZE_GC+1)]
                                if cds_info[presumed_target]['protein_id'] in TARGET_LIST:
                                    params['INC_TARGET_LOADED'] += 1
                                    TARGET_LIST.remove(cds_info[presumed_target]['protein_id'])
                                    cds_info, contig_content = found_target_procedure(presumed_target, cds_info, contig_content)
                    if not TARGET_LIST:
                        break
                ### COM: end of the contig, let the end of the window to treat (target or not, kept or not)
                if TARGET_LIST:
                    window = contig_info[INC_CONTIG_REF]['window']
                    window_length = len(window)
                    contig_content = contig_info[INC_CONTIG_REF]
                    if window_length >= HALF_SIZE_GC:
                        for presumed_target in window[window_length-HALF_SIZE_GC:HALF_SIZE_GC]:
                            if cds_info[presumed_target]['protein_id'] in TARGET_LIST:
                                TARGET_LIST.remove(cds_info[presumed_target]['protein_id'])
                                params['INC_TARGET_LOADED'] += 1
                                cds_info, contig_content = found_target_procedure(presumed_target, cds_info, contig_content)
                        for presumed_target in window[HALF_SIZE_GC:]:
                            if cds_info[presumed_target]['protein_id'] in TARGET_LIST:
                                TARGET_LIST.remove(cds_info[presumed_target]['protein_id'])
                                params['INC_TARGET_LOADED'] += 1
                                cds_info, contig_content = found_target_procedure(presumed_target, cds_info, contig_content)
                            if contig_content['cds_to_keep']:
                                cds_info = remove_useless_cds(window[0], cds_info, contig_content)
                            else:
                                del cds_info[window[0]]
                            del window[0]
                    elif window_length < HALF_SIZE_GC:
                        for presumed_target in window:
                            if cds_info[presumed_target]['protein_id'] in TARGET_LIST:
                                TARGET_LIST.remove(cds_info[presumed_target]['protein_id'])
                                params['INC_TARGET_LOADED'] += 1
                                cds_info, contig_content = found_target_procedure(presumed_target, cds_info, contig_content)
                    logger.debug('proteins in cds_to_keep on this contig: {}'.format(contig_info[INC_CONTIG_REF]['cds_to_keep']))
                    contig_info[INC_CONTIG_REF]['cds_to_keep'] = list(skip_duplicates(contig_info[INC_CONTIG_REF]['cds_to_keep']))
                # COM: completed contig parsing, and still one or more targets have not been found
                if TARGET_LIST:
                    logger.warning('Target(s) has(have) not been found: {}'.format(TARGET_LIST))
            try:
                CONTIG = next(seqRecordParsed)
            except:
                # COM: completed file parsing, and still one or more contigs have not been found
                if CONTIG_LIST:
                    logger.warning('Contig(s) has(have) not been found: {}'.format(CONTIG_LIST))
                break
    return cds_info, contig_info, params

def parse_INSDC_files(d_input, cds_info, contig_info, params):
    ''' every INSDC file in d_input will be parsed
    input: d_input dictionary with filename as key
    output: cds_in and contig_info dictionaries containing information related to contigs
    and targets mentioned in the d_input, and cds belonging to the genomic contexts of mentioned targets
    '''
    logger = logging.getLogger('{}.{}'.format(parse_INSDC_files.__module__, parse_INSDC_files.__name__))
    nbr_of_files = len(d_input)
    for afile in d_input:
        params['INC_FILE'] += 1
        logger.info('Parsing of INSDC file ({}/{}): {}'.format(params['INC_FILE'], nbr_of_files, afile))
        cds_info, contig_info, params = parse_insdc(afile, d_input[afile], cds_info, contig_info, params)
    return cds_info, contig_info, params

def concat_by_dot(alist):
    ''' does the concatenation by a dot
    '''
    return '.'.join(alist)

def write_multiFasta(cds_info, output):
    ''' writes a multiFasta file from the dictionary obtained by the INSDC file
    parser
    '''
    logger = logging.getLogger('{}.{}'.format(write_multiFasta.__module__, write_multiFasta.__name__))
    with open(output, 'w') as fastaFile:
        for cds in cds_info:
            fastaFile.write('>{}\n'.format(cds))
            fastaFile.write(cds_info[cds]['sequence'])
            fastaFile.write('\n')
    if os.path.getsize(output) == 0:
        logger.debug('MultiFasta File is empty.')
        exit(1)
    return 0

def get_lineage(xml):
    ''' extracts the taxonomic lineage from the provided xml file format
    ??? why initialize rank as False and not NA ???
    '''
    logger = logging.getLogger('{}.{}'.format(get_lineage.__module__, get_lineage.__name__))
    lineage_full = {}
    root = ET.fromstring(xml)
    scientificName = root.find('Taxon').find('ScientificName').text
    rank = False
    if root.find('Taxon').find('Rank').text:
        rank = root.find('Taxon').find('Rank').text
    taxId = root.find('Taxon').find('TaxId').text
    collected_taxId = taxId
    lineage_full[scientificName] = {'rank' : rank, 'taxId' : taxId}

    lineage = root.find('Taxon').find('LineageEx')
    taxons = lineage.findall('Taxon')
    for taxon in taxons:
        scientificName = taxon.find('ScientificName').text
        taxId = taxon.find('TaxId').text
        rank = False
        if taxon.find('Rank').text:
            rank = taxon.find('Rank').text
        lineage_full[scientificName] = {'rank' : rank, 'taxId' : taxId}
    return lineage_full, collected_taxId

def get_taxo_from_web(taxonID, tmpDirectoryProcess):
    ''' does web request to get the taxonomic lineage for a taxon ID
    '''
    logger = logging.getLogger('{}.{}'.format(get_taxo_from_web.__module__, get_taxo_from_web.__name__))
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    retry = urllib3.util.Retry(read=5, backoff_factor=2)
    http = urllib3.PoolManager()
    #xml = http.request('GET', 'https://www.ebi.ac.uk/ena/data/view/Taxon:' + str(taxonID) + '&display=xml')
    xml = http.request('GET', 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={}&retmode=xml'.format(str(taxonID)), retries=retry)
    with open('{}/xml_file'.format(tmpDirectoryProcess), 'w') as file:
        file.write(xml.data.decode('utf-8'))
    taxonLineage, collected_taxId = get_lineage(xml.data.decode('utf-8'))
    return taxonLineage, collected_taxId

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
            scientificName = common.global_dict['defaultValue']
            taxId = common.global_dict['defaultValue']
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

def get_taxonLineage(taxonIDs, tmpDirectoryProcess):
    ''' concatenates the various taxonomic lineage dictionaries in a
    super-dictionary called all_lineages
    '''
    logger = logging.getLogger('{}.{}'.format(get_taxonLineage.__module__, get_taxonLineage.__name__))
    logger.info('Getting all taxonomic lineages represented in the Project')
    all_lineages = {}
    for ataxon in taxonIDs:
        res_taxo, collected_taxId = get_taxo_from_web(ataxon, tmpDirectoryProcess)
        if collected_taxId != ataxon:
            logger.info('This species having {} as taxonID, has its taxonID change for {}'.format(ataxon, collected_taxId))
            ataxon = collected_taxId
        desired_taxo = get_desired_lineage(res_taxo)
        new_lineage = store_into_dict(ataxon, desired_taxo)
        all_lineages.update(new_lineage) # no risk of overwriting
    return all_lineages

def taxonomicLineage_runner(contig_info, tmpDirectoryProcess, taxoOut):
    logger = logging.getLogger('{}.{}'.format(taxonomicLineage_runner.__module__, taxonomicLineage_runner.__name__))
    logger.info('Taxonomic lineages research ...')
    taxonIDs = list(set([contig_info[contig]['taxon_ID'] for contig in contig_info if contig_info[contig]['taxon_ID'] != common.global_dict['defaultValue']]))
    taxonomicLineage = get_taxonLineage(taxonIDs, tmpDirectoryProcess)
    logger.info('End of taxonomic lineages research')
    logger.info('Taxonomic lineages information writting ...')
    common.write_pickle(taxonomicLineage, taxoOut)
    common.write_json(taxonomicLineage, '{}/{}'.format(tmpDirectoryProcess, 'taxonomicLineage.json'))
    return 0

def mmseqs_preparation(cds_info, multiFasta, targetsOut):
    ''' write the multifasta file and the targets list
    '''
    write_multiFasta(cds_info, multiFasta)
    targets_storage = [cds for cds in cds_info if cds in cds_info[cds]['target']]
    common.write_pickle(targets_storage, targetsOut)
    common.write_json(targets_storage, '{}/{}/{}'.format(common.global_dict['tmpDirectory'], common.global_dict['boxName']['ClusteringIntoFamilies'], 'targets_list.json'))

def mmseqs_createdb(tmpDirectoryProcess, prefix):
    ''' create database using the mmseqs software
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_createdb.__module__, mmseqs_createdb.__name__))
    suffix = 'faa'
    multiFasta = '{}/{}'.format(tmpDirectoryProcess, concat_by_dot([prefix, suffix]))
    with open('{}/{}'.format(tmpDirectoryProcess, 'mmseqs_createdb.log'), 'w') as file:
        db_creation = subprocess.run(['mmseqs', 'createdb', multiFasta, '{}/{}.{}'.format(tmpDirectoryProcess, prefix, 'DB')], stdout=file, stderr=file, check=True)
        logger.info('createdb - exit code: {}'.format(db_creation.returncode))
    return 0

def mmseqs_clustering(tmpDirectoryProcess, prefix, cov, ident, cov_mode):
    ''' cluster sequences using the mmseqs software
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_clustering.__module__, mmseqs_clustering.__name__))
    suffixDB = 'DB'
    suffixCluster = 'cluster'
    if os.path.isdir('{}/{}'.format(tmpDirectoryProcess, 'MMseqsTMP')):
        shutil.rmtree('{}/{}'.format(tmpDirectoryProcess, 'MMseqsTMP'))
    os.mkdir('{}/{}'.format(tmpDirectoryProcess, 'MMseqsTMP'))
    dataBase = '{}/{}'.format(tmpDirectoryProcess, concat_by_dot([prefix, suffixDB]))
    outputCluster = '{}/{}'.format(tmpDirectoryProcess, concat_by_dot([prefix, suffixCluster]))
    with open('{}/{}'.format(tmpDirectoryProcess, 'mmseqs_clustering.log'), 'w') as file:
        clust_creation = subprocess.run(['mmseqs', 'cluster', dataBase,
                                         outputCluster, '{}/{}'.format(tmpDirectoryProcess, 'MMseqsTMP'),
                                         '--min-seq-id', str(ident),
                                         '--cov-mode', str(cov_mode),
                                         '-c', str(cov),
                                         '--cluster-mode', str(2),
                                         '--kmer-per-seq', str(80),
                                         '--max-seqs', str(300)
                                         ], stdout=file, stderr=file, check=True)
        logger.info('clustering - exit code: {}'.format(clust_creation.returncode))
    return 0

def mmseqs_createTSV(tmpDirectoryProcess, prefix):
    ''' execute the mmseqs command line 'mmseqs createtsv'
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_createTSV.__module__, mmseqs_createTSV.__name__))
    suffixDB = 'DB'
    suffixCluster = 'cluster'
    suffixTSV = 'tsv'
    inputDB = '{}/{}'.format(tmpDirectoryProcess, concat_by_dot([prefix, suffixDB]))
    inputCluster = '{}/{}'.format(tmpDirectoryProcess, concat_by_dot([prefix, suffixCluster]))
    outputTSV = '{}/{}'.format(tmpDirectoryProcess, concat_by_dot([prefix, suffixTSV]))
    with open('{}/{}'.format(tmpDirectoryProcess, 'mmseqs_createtsv.log'), 'w') as file:
        tsv_creation = subprocess.run(['mmseqs', 'createtsv', inputDB,
                                       inputDB, inputCluster, outputTSV
                                       ], stdout=file, stderr=file, check=True)
        logger.info('createTSV - exit code: {}'.format(tsv_creation.returncode))
    return 0

def mmseqs_runner(params, tmpDirectoryProcess):
    ''' runs the mmseqs2 software on the multiFasta file 'ALL.faa'
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_runner.__module__, mmseqs_runner.__name__))
    logger.info('MMseqs2 running ...')
    mmseqs_createdb(tmpDirectoryProcess, params['prefix'])
    mmseqs_clustering(tmpDirectoryProcess, params['prefix'], params['coverage'],
                      params['min_id'], params['cov_mode'])
    mmseqs_createTSV(tmpDirectoryProcess, params['prefix'])
    shutil.rmtree('{}/{}'.format(tmpDirectoryProcess, 'MMseqsTMP/'))
    logger.info('End of MMseqs2 running !')
    return 0

def regroup_families(tsv_file, cds_info):
    ''' create a dictionary to store families obtained by MMseqs2
    '''
    logger = logging.getLogger('{}.{}'.format(regroup_families.__module__, regroup_families.__name__))
    INC_FAMILY = 1
    centroid = None
    lines = common.read_file(tsv_file)
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

def run(INPUT_II, IDENT, COVERAGE):
    ''' main script to run the second box of NetSyn2
    '''
    # Constants
    # common.constantsInitialiszation(args.ProjectName, args.InputFile) # depends on how the function is launched (by hand or via netsyn)
    boxName = common.global_dict['boxName']['ClusteringIntoFamilies']
    tmpDirectoryProcess = '{}/{}'.format(common.global_dict['tmpDirectory'], boxName)
    # Outputs
    multiFasta = common.global_dict['files'][boxName]['faa']
    contigsOut = common.global_dict['files'][boxName]['contigs']
    gcOut = common.global_dict['files'][boxName]['genomicContexts']
    taxoOut = common.global_dict['files'][boxName]['lineage']
    targetsOut = common.global_dict['files'][boxName]['targets']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    logger.info('{} running...'.format(boxName))
    # Process
    if not os.path.isdir(tmpDirectoryProcess):
        os.mkdir(tmpDirectoryProcess)

    params = {
        'PSEUDOGENE': True, # Tells if pseudogenes are included in the analysis
        'MAX_GC': common.global_dict['maxGCSize'],
        'INC_PSEUDO_REF': 0, # counter of pseusogenes
        'INC_CDS_REF': 0,
        'INC_CONTIG_REF': 0,
        'INC_FILE': 0,
        'min_id': IDENT,
        'cov_mode': 1,
        'coverage': COVERAGE,
        }
    params['prefix'] = "MMseqs2_run"

    # MMseq2 files removing (inclure dans netsyn lors nouvelle analyse ?)
    try:
        os.remove('{}/{}.{}'.format(tmpDirectoryProcess, params['prefix'], 'cluster'))
    except:
        pass
    try:
        os.remove('{}/{}.{}'.format(tmpDirectoryProcess, params['prefix'], 'cluster.index'))
    except:
        pass
    try:
        os.remove('{}/{}.{}'.format(tmpDirectoryProcess, params['prefix'], 'tsv'))
    except:
        pass
    try:
        os.remove('{}/{}'.format(tmpDirectoryProcess, 'mmseqs_createtsv.log'))
    except:
        pass

    written_files = os.listdir(tmpDirectoryProcess)
    if ('genomicContexts.pickle' or 'contigs.pickle') not in written_files:
        logger.info('All files have to be written')
        cds_info = {}
        contig_info = {}
        d_input = check_and_get_input(INPUT_II)
        #tester la fonction map() de python pour appliquer une fonction sur une
        #liste
        #usage : map(myFun, myList)
        logger.info('INSDC files parsing ...')
        cds_info, contig_info, params = parse_INSDC_files(d_input, cds_info, contig_info, params)
        logger.info('End of INSDC files parsing !')

        common.write_pickle(contig_info, contigsOut)
        common.write_json(contig_info, '{}/{}'.format(tmpDirectoryProcess, 'contigs.json'))
        common.write_pickle(cds_info, gcOut)
        common.write_json(cds_info, '{}/{}'.format(tmpDirectoryProcess, 'genomicContexts.json'))

        mmseqs_preparation(cds_info, multiFasta, targetsOut)
        exit(1)
        taxonomicLineage_runner(contig_info, tmpDirectoryProcess, taxoOut)

    elif ('MMseqs2_run.faa' or 'targets_list.pickle') not in written_files: # and not ('taxonomicLineage.pickle') ???
        logger.info('Missing the multifasta and targets list files')
        contig_info = common.read_pickle(contigsOut)
        cds_info = common.read_pickle(gcOut)
        mmseqs_preparation(cds_info, multiFasta, targetsOut)
        exit(1)
        taxonomicLineage_runner(contig_info, tmpDirectoryProcess, taxoOut)

    elif 'taxonomicLineage.pickle' not in written_files:
        logger.info('Missing taxonomic lineage file')
        contig_info = common.read_pickle(contigsOut)
        taxonomicLineage_runner(contig_info, tmpDirectoryProcess, taxoOut)

    mmseqs_runner(params, tmpDirectoryProcess)

    try:
        cds_info
    except:
        cds_info = common.read_pickle(gcOut)

    cds_info = regroup_families('{}/{}'.format(tmpDirectoryProcess, concat_by_dot([params["prefix"], 'tsv'])), cds_info)
    common.write_pickle(cds_info, gcOut)
    common.write_json(cds_info, '{}/{}'.format(tmpDirectoryProcess, 'genomicContexts.json'))

    with open('{}/analysed_cds'.format(tmpDirectoryProcess), 'w') as file:
        for inc in cds_info.keys():
            if 'uniprot' in cds_info[inc]:
                uniprot = cds_info[inc]['uniprot']
            else:
                uniprot = common.global_dict['defaultValue']
            file.write('{}\t{}\t{}\n'.format(inc, cds_info[inc]['protein_id'], uniprot))

    logger.info('End of ClusteringIntoFamilies')

if __name__ == '__main__':
    parser = argumentsParser()
    args = parser.parse_args()
    print(args)
    if not os.path.isdir(args.ProjectName):
        os.mkdir(args.ProjectName)
        os.mkdir('{}/TMP'.format(args.ProjectName))
    elif not os.path.isdir('{}/TMP'.format(args.ProjectName)):
        os.mkdir('{}/TMP'.format(args.ProjectName))
    run(args.input, args.Ident, args.MinCoverage)
