#!/usr/bin/env python3

##########
# Import #
##########
import common
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
    ## from JL
    authorized_columns = common.global_dict['inputIIheaders'] + ['taxon_ID']
    mandatory_columns = list(common.global_dict['inputIIheaders'])
    mandatory_columns.remove('UniProt_AC')
    # mandatory_columns = ['protein_AC', 'protein_AC_field', 'nucleic_AC', 'nucleic_File_Format', 'nucleic_File_Name']
    d_rows = parse_tsv(input, authorized_columns, mandatory_columns)
    #read_rows(d_rows)
    ## End from JL
    d_input = create_d_input(d_rows)
    return d_input

def get_from_dbxref(aFeature, dbref):
    ''' retrieve the value from the dbxref list
    '''
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
    taxon_ID = get_uniq_value(aFeature, 'taxon')
    return taxon_ID

def get_orgs_info(aFeature, taxon_ID, numFile):
    ''' get the relied information to a contig as organism, strain, size
    and taxon ID
    input: source feature from an INSDC file
    output: updated contig dictionary
    '''
    if taxon_ID == common.global_dict['defaultValue']:
        taxon_ID = get_taxonID(aFeature)
    org_content = {
        'id': numFile,
        'name': get_uniq_value(aFeature, 'organism'),
        'strain': get_uniq_value(aFeature, 'strain'),
        'taxon_ID': taxon_ID
        # 'size': [aFeature.location.start.real+1,
        #          aFeature.location.end.real
        #          ]
        }
    return org_content

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

# def det_frame(direction, startCDS, endGenome):
#     ''' determine the frame of the CDS
#     -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
#     -*-*-*- needs to be rethought because of cases where there is a   -*-*-*-
#     -*-*-*- join(##..##, ##..###)                                     -*-*-*-
#     -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
#     '''
#     if direction == '1':
#         window = ((int(startCDS)-1)%3)+1
#         strand = ''.join(['+', str(window)])
#     elif direction == '-1':
#         window = ((int(endGenome)-int(startCDS))%3)+1
#         strand = ''.join(['-', str(window)])
#     return strand

def get_pseudo_info(aFeature, sequences, window, params):
    ''' add information of the pseudogene to prots_info dictionary
    input: CDS feature having a '/pseudo' qualifier
    output: updated prots_info dictionary, window list of the contig_content
    modified to contain afeature
    '''
    params['INC_CDS_REF'] += 1
    INC_CDS_REF = params['INC_CDS_REF']

    pseudo_id = get_pseudo_id(params)
    begin, end = aFeature.location.start.real+1, aFeature.location.end.real

    feature_info= {
        'id': INC_CDS_REF,
        'protein_id': pseudo_id,
        'begin': begin,
        'end': end,
        'strand': str(aFeature.location.strand),
        'product': aFeature.qualifiers.get('product') if aFeature.qualifiers.get('product') else common.global_dict['defaultValue'],
        'targets': []
        }

    sequences[INC_CDS_REF] = get_uniq_value(aFeature, 'translation')
    window.append(INC_CDS_REF)
    return feature_info, sequences, window, params

def get_prot_info(aFeature, sequences, window, proteinField, params):
    ''' add information of the protein to prots_info dictionary
    input: CDS feature of an INSDC file
    output: updated prots_info dictionary, window list of the contig_content
    modified to contain afeature
    '''
    logger = logging.getLogger('{}.{}'.format(get_prot_info.__module__, get_prot_info.__name__))
    params['INC_CDS_REF'] += 1
    INC_CDS_REF = params['INC_CDS_REF']

    ident = get_required_value(get_uniq_value, aFeature, proteinField)
    if ident == 1:
        ident = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(15))
        logger.info('Protein {} has no field protein_id. A random name has been created: {}'.format(INC_CDS_REF, ident))
    begin, end = aFeature.location.start.real+1, aFeature.location.end.real
    ec_nums = (aFeature.qualifiers.get('EC_number') if aFeature.qualifiers.get('EC_number') else common.global_dict['defaultValue'])

    feature_info = {
        'id': INC_CDS_REF,
        'protein_id': ident,
        'begin': begin,
        'end': end,
        'strand': str(aFeature.location.strand),
        'product': aFeature.qualifiers.get('product') if aFeature.qualifiers.get('product') else common.global_dict['defaultValue'],
        'ec_number': ec_nums,
        'uniprot': get_from_dbxref(aFeature, 'UniProt'),
        'targets': []
        #'gene_name'
        #'locus_tag'
        }

    sequences[INC_CDS_REF] = get_required_value(get_uniq_value, aFeature, 'translation'),
    window.append(INC_CDS_REF)
    return feature_info, sequences, window, params

def set_target_to_gc(ref_target, gcList, prots_info):
    ''' add the target reference to members of its genomic context
    input: target reference, list of target's genomic context
    output: updated prots_info dictionary on genomic context cds key
    '''
    for prot_dict in prots_info[-common.global_dict['maxGCSize']:]:
        prot_idx = prots_info.index(prot_dict)
        if prots_info[prot_idx]['id'] in gcList:
        #if ref_target not in prots_info[gc_member]['target']:
            prots_info[prot_idx]['targets'].append(ref_target)
    return prots_info

def found_target_procedure(target, prots_info, targets_info, cds_to_keep, window):
    ''' add the genomic context to the 'cds_to_keep' list, copy the 'window'
    to the 'context'
    input: target reference
    output: updated prots_info on target's context and updated contig_content on
    cds_to_keep
    '''
    logger = logging.getLogger('{}.{}'.format(found_target_procedure.__module__, found_target_procedure.__name__))
    cds_to_keep.extend(window)

    targets_info.append({
            'id': target,
            'context': window.copy()
            #'organism_id': ,
            })
    prots_info = set_target_to_gc(target, window, prots_info)
    #logger.debug('target ({}/{})\t- contig_content ({}) -\twindow: {}'.format(target, prots_info[target]['protein_id'], contig_content['contig'], contig_content['window']))
    #logger.debug('target ({}/{})\t- contig_content ({}) -\tcds_to_keep: {}'.format(target, prots_info[target]['protein_id'], contig_content['contig'], contig_content['cds_to_keep']))
    return prots_info, cds_to_keep

def is_useless_cds(cds, cds_to_keep):
    ''' remove the cds if is not in cds_to_keep
    input: cds to test
    output: potential modified prots_info dictionary where cds entry is removed
    '''
    return cds not in cds_to_keep

def parse_insdc(afile, d_infile, prots_info, targets_info, orgs_info, sequences, params):
    ''' get information for every contig mentionned in the input
    and get information for every cds that is contained in a (MAX_GC sized-)
    window centered on a target
    input: INSDC file, file input information in d_infile dictionary
    output: prots_info dictionary containing the cds belonging to the genomic context
    of the targets and targets itself, and
    orgs_info dictionary containing the explored contig information
    (see get_orgs_info function), and
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
                #orgs_info[INC_CONTIG_REF] = {'contig': contig_name}
                TARGET_LIST = d_infile[contig_name]['target_list']
                cds_to_keep = []
                window = []
                logger.debug('Target list for {} contig: {}'.format(contig_name, TARGET_LIST))

                # COM: beginning of the parsing
                for aFeature in CONTIG.features:
                    if aFeature.type == 'source': ####### il faut encore tester si on a déjà récupéré les infos dans le cas d'un fichier INSDC multicontig
                        orgs_info.append(
                            get_orgs_info(
                                aFeature,
                                d_infile[contig_name]['taxon_ID']
                                if 'taxon_ID' in d_infile[contig_name].keys()
                                else common.global_dict['defaultValue'],
                                params['INC_FILE']
                                )
                            )
                        print(orgs_info)
                    elif aFeature.type == 'CDS':
                        newCdsAdded = False
                        if is_pseudogene(aFeature):
                            if params['PSEUDOGENE']:
                                feature_info, sequences, window, params = get_pseudo_info(aFeature, sequences, window, params)
                                prots_info.append(feature_info)
                                newCdsAdded = True
                        else:
                            feature_info, sequences, window, params = get_prot_info(aFeature, sequences, window, fieldProteinID, params)
                            prots_info.append(feature_info)
                            newCdsAdded = True

                        if newCdsAdded:
                            window_length = len(window)
                            if window_length == MAX_GC:
                                presumed_target = window[HALF_SIZE_GC]
                                presumed_kicked_out_cds = window[0]

                                if prots_info[-(HALF_SIZE_GC+1)]['protein_id'] in TARGET_LIST:
                                    params['INC_TARGET_LOADED'] += 1
                                    TARGET_LIST.remove(prots_info[-(HALF_SIZE_GC+1)]['protein_id'])
                                    prots_info, cds_to_keep = found_target_procedure(presumed_target, prots_info, targets_info, cds_to_keep, window)
                                if cds_to_keep:
                                    if is_useless_cds(presumed_kicked_out_cds, cds_to_keep):
                                        del prots_info[-MAX_GC]
                                else:
                                    del prots_info[-MAX_GC]
                                del window[0]
                            elif window_length > HALF_SIZE_GC:
                                presumed_target = window[-(HALF_SIZE_GC+1)]
                                
                                if prots_info[-(HALF_SIZE_GC+1)]['protein_id'] in TARGET_LIST:
                                    params['INC_TARGET_LOADED'] += 1
                                    TARGET_LIST.remove(prots_info[-(HALF_SIZE_GC+1)]['protein_id'])
                                    prots_info, cds_to_keep = found_target_procedure(presumed_target, prots_info, targets_info, cds_to_keep, window)

                    if not TARGET_LIST:
                        break
                ### COM: end of the contig, let the end of the window to treat (target or not, kept or not)
                if TARGET_LIST:
                    window_length = len(window)
                    if window_length >= HALF_SIZE_GC:
                        while window_length > HALF_SIZE_GC:
                            if prots_info[-(window_length-HALF_SIZE_GC)]['protein_id'] in TARGET_LIST:
                                TARGET_LIST.remove(prots_info[-(window_length-HALF_SIZE_GC)]['protein_id'])
                                params['INC_TARGET_LOADED'] += 1
                                prots_info, cds_to_keep = found_target_procedure(presumed_target, prots_info, targets_info, cds_to_keep, window)
                            if cds_to_keep:
                                if is_useless_cds(prots_info[-window_length]['id'], cds_to_keep):
                                    del prots_info[-window_length]
                            else:
                                del prots_info[-window_length]
                            del window[0]
                            window_length = len(window)
                    # elif window_length < HALF_SIZE_GC:
                    #     for presumed_target in window:
                    #         if prots_info[presumed_target]['protein_id'] in TARGET_LIST:
                    #             TARGET_LIST.remove(prots_info[presumed_target]['protein_id'])
                    #             params['INC_TARGET_LOADED'] += 1
                    #             prots_info, contig_content = found_target_procedure(presumed_target, prots_info, contig_content)

                    #logger.debug('proteins in cds_to_keep on this contig: {}'.format(orgs_info[INC_CONTIG_REF]['cds_to_keep']))
                    #orgs_info[INC_CONTIG_REF]['cds_to_keep'] = list(skip_duplicates(orgs_info[INC_CONTIG_REF]['cds_to_keep']))
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
    return prots_info, targets_info, orgs_info, sequences, params

def parse_INSDC_files(d_input, prots_info, targets_info, orgs_info, params):
    ''' every INSDC file in d_input will be parsed
    input: d_input dictionary with filename as key
    output: cds_in and orgs_info dictionaries containing information related to contigs
    and targets mentioned in the d_input, and cds belonging to the genomic contexts of mentioned targets
    '''
    logger = logging.getLogger('{}.{}'.format(parse_INSDC_files.__module__, parse_INSDC_files.__name__))
    nbr_of_files = len(d_input)
    sequences = {}
    for afile in d_input:
        params['INC_FILE'] += 1
        logger.info('Parsing of INSDC file ({}/{}): {}'.format(params['INC_FILE'], nbr_of_files, afile))
        prots_info, targets_info, orgs_info, sequences, params = parse_insdc(afile, d_input[afile], prots_info, targets_info, orgs_info, sequences, params)
    exit(1)
    return prots_info, targets_info, orgs_info, params

def concat_by_dot(alist):
    ''' does the concatenation by a dot
    '''
    return '.'.join(alist)

def write_multiFasta(prots_info, output):
    ''' writes a multiFasta file from the dictionary obtained by the INSDC file
    parser
    '''
    logger = logging.getLogger('{}.{}'.format(write_multiFasta.__module__, write_multiFasta.__name__))
    with open(output, 'w') as fastaFile:
        for cds in prots_info:
            fastaFile.write('>{}\n'.format(cds))
            fastaFile.write(prots_info[cds]['sequence'])
            fastaFile.write('\n')
    if os.path.getsize(output) == 0:
        logger.debug('MultiFasta File is empty.')
        exit(1)
    return 0


def get_lineage(xml, desiredTaxonIDs):
    logger = logging.getLogger('{}.{}'.format(get_lineage.__module__, get_lineage.__name__))
    ''' extracts the taxonomic lineage from the provided xml file format
    ??? why initialize rank as False and not NA ???
    '''
    allLineages = []
    root = ET.fromstring(xml)
    if root.findall('Taxon'):
        for taxon in root.findall('Taxon'):
            fullLineage = {}
            collectedTaxId = taxon.find('TaxId').text
            scientificName = taxon.find('ScientificName').text
            rank = False
            if taxon.find('Rank').text:
                rank = root.find('Taxon').find('Rank').text
            fullLineage[scientificName] = {
                                    'rank' : rank,
                                    'taxId' : collectedTaxId
                                    }
            for subTaxon in taxon.find('LineageEx').findall('Taxon'):
                taxId = subTaxon.find('TaxId').text
                scientificName = subTaxon.find('ScientificName').text
                rank = False
                if subTaxon.find('Rank').text:
                    rank = subTaxon.find('Rank').text
                fullLineage[scientificName] = {
                                    'rank' : rank,
                                    'taxId' : taxId
                                    }
            oldTaxonIDs = []
            if taxon.find('AkaTaxIds'):
                for taxon in taxon.find('AkaTaxIds').findall('TaxId'):
                    oldTaxonIDs.append(taxon.text)
            allLineages.append([collectedTaxId, fullLineage, oldTaxonIDs])
    indexTaxID = 0
    # indexLineage = 1
    indexOldTaxID = 2
    if len(desiredTaxonIDs) == len(allLineages):
        # All lineage recovered
        for index, desiredTaxon in enumerate(desiredTaxonIDs):
            collectedTaxId = allLineages[index][indexTaxID]
            if collectedTaxId == desiredTaxon:
                logger.info('Taxonomic lineage recovered for the taxonID {}.'.format(desiredTaxon))
            elif desiredTaxon in allLineages[index][indexOldTaxID]:
                logger.info('The species having {} as taxonID, has its taxonID change for {}'.format(desiredTaxon, collectedTaxId))
                logger.warning('§§§§§ IL FAUT METTRE A JOUR LE TAXID DES ORGANISMES CORRESPONDANTS DANS CONTIGS §§§§§')
    elif len(allLineages) > 0:
        # Lineage missing
        indexLineage = 0
        indexDesiredTaxon = 0
        while 1:
            collectedTaxId = allLineages[indexLineage][indexTaxID]
            oldTaxonIDs = allLineages[indexLineage][indexOldTaxID]
            desiredTaxon = desiredTaxonIDs[indexDesiredTaxon]
            if collectedTaxId == desiredTaxon:
                logger.info('Taxonomic lineage recovered for the taxonID {}.'.format(desiredTaxon))
                indexLineage += 1
                indexDesiredTaxon += 1
            elif desiredTaxon in allLineages[indexLineage][indexOldTaxID]:
                logger.info('The species having {} as taxonID, has its taxonID change for {}'.format(desiredTaxon, collectedTaxId))
                logger.warning('§§§§§ IL FAUT METTRE A JOUR LE TAXID DES ORGANISMES CORRESPONDANTS DANS CONTIGS §§§§§')
                indexLineage += 1
                indexDesiredTaxon += 1
            else:
                logger.warning('NCBI doesn\'t reconise the taxonID: {}.'.format(desiredTaxon))
                allLineages.append([desiredTaxon, {}, []])
                indexDesiredTaxon += 1
            if indexDesiredTaxon == len(desiredTaxonIDs):
                break
            elif indexLineage == len(allLineages):
                for desiredTaxon in desiredTaxonIDs[indexDesiredTaxon:]:
                    logger.warning('NCBI doesn\'t reconise the taxonID: {}.'.format(desiredTaxon))
                    allLineages.append([desiredTaxon, {}, []])
                break
    else:
        for desiredTaxon in desiredTaxonIDs:
            logger.warning('NCBI doesn\'t reconise the taxonID: {}.'.format(desiredTaxon))
    return allLineages

def get_taxo_from_web(taxonIDs, tmpDirectoryProcess):
    ''' does web request to get the taxonomic lineage for a taxon ID
    '''
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    http = urllib3.PoolManager()
    ids = ','.join(taxonIDs)
    xml = common.httpRequest(http,'GET', 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={}&retmode=xml'.format(ids))
    with open('{}/xml_file'.format(tmpDirectoryProcess), 'w') as file:
        file.write(xml.data.decode('utf-8'))
    return get_lineage(xml.data.decode('utf-8'), taxonIDs)

def get_desired_lineage(lineages):
    ''' takes the ranks of interest through the taxonomic lineage
    '''
    allLineages = {}
    for collectedTaxId, fullLineage, oldTaxonIDs in lineages:
        desired_ranks = dict(common.global_dict['desired_ranks_lneage'])
        desired_lineage = []
        for scientificName in fullLineage:
            if fullLineage[scientificName]['rank'] in desired_ranks:
                rank = fullLineage[scientificName]['rank']
                level = desired_ranks[rank]
                taxId = fullLineage[scientificName]['taxId']
                desired_lineage.append([level, rank, scientificName, taxId])
                del desired_ranks[rank]
        if desired_ranks != {}:
            for rank in desired_ranks:
                level = desired_ranks[rank]
                scientificName = 'NA'
                taxId = 'NA'
                desired_lineage.append([level, rank, scientificName, taxId])
        allLineages.update(store_into_dict(collectedTaxId, desired_lineage))
    return allLineages

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

def get_taxonomicLineage(taxonIDs, tmpDirectoryProcess):
    ''' concatenates the various taxonomic lineage dictionaries in a
    super-dictionary called all_lineages
    '''
    logger = logging.getLogger('{}.{}'.format(get_taxonomicLineage.__module__, get_taxonomicLineage.__name__))
    logger.info('Getting all taxonomic lineages represented in the Project')
    all_lineages = {}
    batchesSize = 100
    # taxonIDs = ['1884790', '2', '3', '304371', '3','4', '76869','3','4','3','4','3','4','3','4'] # To test
    while taxonIDs:
        lineages = get_taxo_from_web(taxonIDs[:batchesSize], tmpDirectoryProcess)
        all_lineages.update(get_desired_lineage(lineages))
        del taxonIDs[:batchesSize]
    return all_lineages

def taxonomicLineage_runner(orgs_info, tmpDirectoryProcess, taxoOut):
    logger = logging.getLogger('{}.{}'.format(taxonomicLineage_runner.__module__, taxonomicLineage_runner.__name__))
    logger.info('Taxonomic lineages research ...')
    taxonIDs = list(set([orgs_info[contig]['taxon_ID'] for contig in orgs_info if orgs_info[contig]['taxon_ID'] != common.global_dict['defaultValue']]))
    taxonomicLineage = get_taxonomicLineage(taxonIDs, tmpDirectoryProcess)
    logger.info('End of taxonomic lineages research')
    logger.info('Taxonomic lineages information writting ...')
    common.write_pickle(taxonomicLineage, taxoOut)
    common.write_json(taxonomicLineage, '{}/{}'.format(tmpDirectoryProcess, 'taxonomicLineage.json'))
    return 0

def mmseqs_preparation(prots_info, multiFasta, targetsOut, tmpDirectoryProcess):
    ''' write the multifasta file and the targets list
    '''
    write_multiFasta(prots_info, multiFasta)
    targets_storage = [cds for cds in prots_info if cds in prots_info[cds]['target']]
    common.write_pickle(targets_storage, targetsOut)
    common.write_json(targets_storage, '{}/{}'.format(tmpDirectoryProcess, 'targets_list.json'))

def mmseqs_createdb(tmpDirectoryProcess, multiFasta, prefix):
    ''' create database using the mmseqs software
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_createdb.__module__, mmseqs_createdb.__name__))
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

def mmseqs_runner(params, tmpDirectoryProcess, multiFasta):
    ''' runs the mmseqs2 software on the multiFasta file 'ALL.faa'
    '''
    logger = logging.getLogger('{}.{}'.format(mmseqs_runner.__module__, mmseqs_runner.__name__))
    logger.info('MMseqs2 running ...')
    mmseqs_createdb(tmpDirectoryProcess, multiFasta, params['prefix'])
    mmseqs_clustering(tmpDirectoryProcess, params['prefix'], params['coverage'],
                      params['min_id'], params['cov_mode'])
    mmseqs_createTSV(tmpDirectoryProcess, params['prefix'])
    shutil.rmtree('{}/{}'.format(tmpDirectoryProcess, 'MMseqsTMP/'))
    logger.info('End of MMseqs2 running !')
    return 0

def regroup_families(tsv_file, prots_info):
    ''' create a dictionary to store families obtained by MMseqs2
    '''
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
        prots_info[cds]['similarityFamily'] = INC_FAMILY
    return prots_info

def run(INPUT_II, insdcDirectory, IDENT, COVERAGE):
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
        'PSEUDOGENE': False, # Tells if pseudogenes are included in the analysis
        'MAX_GC': common.global_dict['maxGCSize'],
        'INC_PSEUDO_REF': 0, # counter of pseusogenes
        'INC_CDS_REF': 0,
        'INC_CONTIG_REF': 0,
        'INC_FILE': 0,
        'min_id': IDENT,
        'cov_mode': 1,
        'coverage': COVERAGE,
        }
    params['prefix'] = multiFasta.split('.')[0]

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
        prots_info = []
        targets_info = []
        orgs_info = []
        d_input = check_and_get_input(INPUT_II)
        #tester la fonction map() de python pour appliquer une fonction sur une
        #liste
        #usage : map(myFun, myList)
        logger.info('INSDC files parsing ...')
        prots_info, targets_info, orgs_info, params = parse_INSDC_files(d_input, prots_info, targets_info, orgs_info, params)
        logger.info('End of INSDC files parsing !')

        common.write_pickle(orgs_info, contigsOut)
        common.write_json(orgs_info, '{}/{}'.format(tmpDirectoryProcess, 'contigs.json'))
        common.write_pickle(prots_info, gcOut)
        common.write_json(prots_info, '{}/{}'.format(tmpDirectoryProcess, 'genomicContexts.json'))

        mmseqs_preparation(prots_info, multiFasta, targetsOut, tmpDirectoryProcess)
        taxonomicLineage_runner(orgs_info, tmpDirectoryProcess, taxoOut)

    elif ('MMseqs2_run.faa' or 'targets_list.pickle') not in written_files: # and not ('taxonomicLineage.pickle') ???
        logger.info('Missing the multifasta and targets list files')
        orgs_info = common.read_pickle(contigsOut)
        prots_info = common.read_pickle(gcOut)
        mmseqs_preparation(prots_info, multiFasta, targetsOut,tmpDirectoryProcess)
        taxonomicLineage_runner(orgs_info, tmpDirectoryProcess, taxoOut)

    elif 'taxonomicLineage.pickle' not in written_files:
        logger.info('Missing taxonomic lineage file')
        orgs_info = common.read_pickle(contigsOut)
        taxonomicLineage_runner(orgs_info, tmpDirectoryProcess, taxoOut)

    mmseqs_runner(params, tmpDirectoryProcess, multiFasta)

    try:
        prots_info
    except:
        prots_info = common.read_pickle(gcOut)

    prots_info = regroup_families('{}/{}'.format(tmpDirectoryProcess, concat_by_dot([params["prefix"], 'tsv'])), prots_info)
    common.write_pickle(prots_info, gcOut)
    common.write_json(prots_info, '{}/{}'.format(tmpDirectoryProcess, 'genomicContexts.json'))

    # with open('{}/analysed_cds'.format(tmpDirectoryProcess), 'w') as file:
    #     for inc in prots_info.keys():
    #         if 'uniprot' in prots_info[inc]:
    #             uniprot = prots_info[inc]['uniprot']
    #         else:
    #             uniprot = common.global_dict['defaultValue']
    #         file.write('{}\t{}\t{}\n'.format(inc, prots_info[inc]['protein_id'], uniprot))

    logger.info('End of ClusteringIntoFamilies')

def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(usage='''ClusteringIntoFamilies.py -i <inputFileName> -o <OutputName> -insdc <insdcDirectoryName> -id <ident> -mc <MinimuCoverage>''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-i', '--input', type=str,
                        required=True, help='File of corresponding.')
    group1.add_argument('-o', '--OutputName', type=str,
                        required=True, help='Output name files.')
    group1.add_argument('-insdc', '--insdcDirectory', type = str,
                        help = 'Directory containing the INSDC files.')
    group1.add_argument('-id', '--Ident', type=float,
                        default=0.3, help='Sequence identity.\nDefault value: 0.3.')
    group1.add_argument('-mc', '--MinCoverage', type=float,
                        default=0.8, help='Minimal coverage allowed.\nDefault value: 0.8.')

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
    #########################
    # Dependancies checking #
    #########################
    common.dependanciesChecking()
    #############
    # Constants #
    #############
    common.global_dict['tmpDirectory'] = '.'
    boxName = common.global_dict['boxName']['ClusteringIntoFamilies']
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('faa', '{}.faa'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('contigs', '{}_contigs.pickle'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('genomicContexts', '{}_genomicsContexts.pickle'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('lineage', '{}_lineage.pickle'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('targets', '{}_targets.pickle'.format(args.OutputName))
    #######
    # Run #
    #######
    run(args.input, args.insdcDirectory ,args.Ident, args.MinCoverage)
