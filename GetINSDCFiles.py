#!/usr/bin/env python3

##########
# Import #
##########
import os
import shutil
import logging
import urllib3
import gzip
import re
import common

#############
# Functions #
#############
def resultsFormat(res, dico):
    '''
    Formating ID match and update dico.
    '''
    logger = logging.getLogger('{}.{}'.format(resultsFormat.__module__, resultsFormat.__name__))
    sepColumn = '\t'
    sepField = ';'
    isHeader = True
    for bLine in res.data.splitlines():
        if isHeader:
            headers = bLine.decode('utf-8').split(sepColumn)
            if not headers[0] == 'Entry':
                logger.critical('"id" column must be in first position')
                exit(2)
            isHeader = False
        else:
            for index, column in enumerate(bLine.decode('utf-8').split(sepColumn)):
                if index == 0:
                    entry = column
                    dico[entry] = {}
                elif not column:
                    logger.warning('{}: This entry is obsolete or not found. Please check in www.UniProt.org.'.format(entry))
                    del dico[entry]
                    break
                else:
                    dico[entry][headers[index]] = column.strip(sepField).split(sepField)
    return dico

def getENAidMatchingToUniProtid(uniprotAccessions, batchesSize, PoolManager):
    '''
    Allows the correspondence between a UniProt accession and nuclotide accessions.
    Batch splitting to lighten the request.
    '''
    logger = logging.getLogger('{}.{}'.format(getENAidMatchingToUniProtid.__module__, getENAidMatchingToUniProtid.__name__))
    crossReference = {}
    nbTotalEntries = len(uniprotAccessions)
    nbEntriesProcessed = 0
    logger.info('Beginning of the correspondence between the UniProt and ENA identifiers...')
    while uniprotAccessions:
        accessions = '+OR+id:'.join(uniprotAccessions[:batchesSize])
        res = common.httpRequest(PoolManager,'GET','https://www.uniprot.org/uniprot/?query=id:{}&columns=id,database(EMBL),database(EMBL_CDS)&format=tab'.format(accessions))
        crossReference = resultsFormat(res, crossReference)
        nbEntriesProcessed += len(uniprotAccessions[:batchesSize])
        del uniprotAccessions[:batchesSize]
        logger.info('Correspondences computed: {}/{}'.format(nbEntriesProcessed,nbTotalEntries))
    return crossReference

def getNucleicFileName(nucleicAccession):
    '''
    Get nucleic file name from nucleic accession.
    '''
    nucleicAcc = re.match(r'^(?P<NucleicFileName>[A-Z]{4,6}[0-9]{2})[0-9]{6,}$', nucleicAccession)
    if nucleicAcc:
        return nucleicAcc.group("NucleicFileName")
    return nucleicAccession

def getEMBLfromENA(nucleicAccession, nucleicFilePath, PoolManager):
    '''
    Download EMBL file from ENA website.
    '''
    logger = logging.getLogger('{}.{}'.format(getEMBLfromENA.__module__, getEMBLfromENA.__name__))
    res = common.httpRequest(PoolManager, 'GET', 'https://www.ebi.ac.uk/ena/data/view/{}&display=text&set=true'.format(nucleicAccession))
    contentType = res.info()['Content-Type']
    if contentType == 'text/plain;charset=UTF-8' and res.data.decode('utf-8') == 'Entry: {} display type is either not supported or entry is not found.\n'.format(nucleicAccession):
        logger.error(res.data.decode('utf-8'))
        exit(1)
    with open(nucleicFilePath, 'w') as file:
        if contentType == 'text/plain;charset=UTF-8':
            file.write(res.data.decode('utf-8'))
        elif contentType == 'application/x-gzip':
            file.write(gzip.decompress(res.data).decode('utf-8'))
        else:
          logger.critical('Unsupported content type ({}).'.format(contentType))
          exit(1)
    logger.info('{} downloaded.'.format(nucleicFilePath))

def run(InputName):
    '''
    Get INSDC files porocessing.
    '''
    # Constants
    boxName = common.global_dict['boxName']['GetINSDCFiles']
    dataDirectoryProcess = '{}/{}'.format(common.global_dict['dataDirectory'], boxName)
    outputName = common.global_dict['files'][boxName]['inputClusteringStep']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    print('')
    logger.info('{} running...'.format(boxName))
    # Process
    if not os.path.isdir(dataDirectoryProcess):
        os.mkdir(dataDirectoryProcess)
    if os.path.isfile(outputName):
        os.remove(outputName)
    header, accessions = common.parseInputI(InputName)
    if header == common.global_dict['inputIheader']:
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        http = urllib3.PoolManager()
        crossReference = getENAidMatchingToUniProtid(accessions, 500, http)
        outputContent = []
        logger.info('Beginning of EMBL file downloading...')
        for entry in crossReference:
            maxAssemblyLength = 0
            for index, nucleicAccession in enumerate(crossReference[entry]['Cross-reference (EMBL)']):
                filesExtension = common.global_dict['filesExtension']
                nucleicFilePath = '{}/{}.{}'.format(dataDirectoryProcess, getNucleicFileName(nucleicAccession), filesExtension)
                if not os.path.isfile(nucleicFilePath):
                    getEMBLfromENA(nucleicAccession, nucleicFilePath, http)
                else:
                    logger.info('{} already existing.'.format(nucleicFilePath))
                with open(nucleicFilePath, 'r') as file:
                    m = re.match(r'^ID .*; (?P<length>[0-9]+) BP\.', file.read())
                    if m:
                        assemblyLength = int(m.group('length'))
                    else:
                        logger.error('{}: improper file format.'.format(nucleicFilePath))
                        exit(1)
                if assemblyLength > maxAssemblyLength:
                    maxAssemblyLength = assemblyLength
                    toAppend = [
                        entry,
                        crossReference[entry]['Cross-reference (embl)'][index],
                        'protein_id',
                        nucleicAccession,
                        'embl',
                        nucleicFilePath
                    ]
            outputContent.append(toAppend)
        with open(outputName, 'w') as file:
            file.write('{}\n'.format('\t'.join(common.global_dict['inputIIheaders'])))
            for line in ['\t'.join(values) for values in outputContent]:
                file.write('{}\n'.format(line))
        logger.info('{} generated.'.format(outputName))
        logger.info('{} completed!'.format(boxName))
    else:
        logger.error('Input header unrecognized.')
        exit(1)

def argumentsParser():
    '''
    Arguments parsing.
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                 usage = '''GetINSDCFiles.py -u <UniProtAC.list> -o <OutputName>''', ######################################################################
                                 formatter_class = argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-u', '--UniProtACList', type = str,
                        required = True, help = 'Protein accession list.')
    group1.add_argument('-o', '--OutputName', type = str,
                        required = True, help = 'Name of corresponding file.')

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
    boxName = common.global_dict['boxName']['GetINSDCFiles']
    common.global_dict.setdefault('files', {}).setdefault(boxName,{}).setdefault('inputClusteringStep', args.OutputName)
    #######
    # Run #
    #######
    run(args.UniProtACList)
