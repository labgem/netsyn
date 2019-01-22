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
def parseInputI(filename): # Fonction a deplacer dans tools ??
    '''
    Input file parssing.
    '''
    logger = logging.getLogger('{}.{}'.format(parseInputI.__module__, parseInputI.__name__))
    with open(filename, 'r') as file:
        firstLine = True
        accessions = []
        seps = [' ','\t',',',';']
        for line in file:
            for sep in seps:
                if len(line.split(sep)) > 1:
                    logger.error('Input invalidate: unauthorized character {}'.format(seps))
            if firstLine:
                header = line.rstrip()
                firstLine = False
            else:
                accessions.append(line.rstrip())
    return header, accessions

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
                logger.critical('"id" column should be in first position.')
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
    while uniprotAccessions:
        accessions = '+OR+id:'.join(uniprotAccessions[:batchesSize])
        resStatus = 0
        remainingTry = 5
        while not resStatus == 200 and remainingTry > 0:
            try:
                logger.info('Matching between UniProt and ENA. Connection to https://www.uniprot.org/, remaining try: {}.'.format(remainingTry))
                res = PoolManager.request('GET' ,
                    'https://www.uniprot.org/uniprot/?query=id:{}&columns=id,database(EMBL),database(EMBL_CDS)&format=tab'.format(accessions))
                resStatus = res.status
                remainingTry -= 1
            except:
                logger.error('OUPS')
                exit(1)
        if not res.status == 200:
            logger.error('HTTP error {}! Matching between UniProt and ENA failed.'.format(res.status))
            exit(1)
        else:
            logger.info('Connection to https://www.uniprot.org/, success!')
        crossReference = resultsFormat(res, crossReference)
        del uniprotAccessions[:batchesSize]
    return crossReference

def getNucleicFialeName(nucleicAccession):
    nucleicAcc = re.match(r'(?P<NucleicFileName>.{6})[0-9]{6}[0-9]*$', nucleicAccession)
    #r'(?P<NucleicFileName>[A-B]{4,6}[0-1]{2})[0-9]{6,}$' ####### revenir dessus plus tard
    if nucleicAcc:
        return nucleicAcc.group("NucleicFileName")
    return nucleicAccession

def getEMBLfromENA(nucleicAccession, nucleicFilePath, PoolManager):
    '''
    Download EMBL file from ENA website.
    '''
    logger = logging.getLogger('{}.{}'.format(getEMBLfromENA.__module__, getEMBLfromENA.__name__))
    resStatus = 0
    remainingTry = 5
    while not resStatus == 200 and remainingTry > 0:
        try:
            logger.info('{} dowloading. Connection to https://www.ebi.ac.uk/ena, remaining try: {}.'.format(nucleicFilePath, remainingTry))
            res = PoolManager.request('GET' , 'https://www.ebi.ac.uk/ena/data/view/{}&display=text&set=true'.format(nucleicAccession))
            resStatus = res.status
            remainingTry -= 1
        except:
            logger.error('OUPS')
            exit(1)
    contentType = res.info()['Content-Type']
    if not res.status == 200:
        logger.error('HTTP error {}! {}.embl couldn\'t be downloaded !'.format(res.status, nucleicFilePath))
        exit(1)
    elif contentType == 'text/plain;charset=UTF-8' and res.data.decode('utf-8') == 'Entry: {} display type is either not supported or entry is not found.\n'.format(nucleicAccession):
        logger.error(res.data.decode('utf-8'))
        exit(1)
    else:
        logger.info('Connection to https://www.uniprot.org/, success!')
    with open(nucleicFilePath, 'w') as file:
        if contentType == 'text/plain;charset=UTF-8':
            file.write(res.data.decode('utf-8'))
        elif contentType == 'application/x-gzip':
            file.write(gzip.decompress(res.data).decode('utf-8'))
        else:
          logger.critical('Unsupported content type ({}).'.format(contentType))
          exit(1)

def run(InputName):
    '''
    Get INSDC files porocessing.
    '''
    # Constants
    boxName = common.global_dict['boxName']['GetINSDCFiles']
    tmpDirectoryProcess = '{}/{}'.format(common.global_dict['tmpDirectory'], boxName)
    outputName = common.global_dict['files'][boxName]['inputClusteringStep']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    logger.info('{} running...'.format(boxName))
    # Process
    if not os.path.isdir(tmpDirectoryProcess):
        os.mkdir(tmpDirectoryProcess)
    if os.path.isfile(outputName):
        os.remove(outputName)
    header, accessions = parseInputI(InputName)
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    http = urllib3.PoolManager()
    if header == 'UniProtAC':
        crossReference = getENAidMatchingToUniProtid(accessions, 500, http)
        outputContent = []
        for entry in crossReference:
            maxAssemblyLength = 0
            for index, nucleicAccession in enumerate(crossReference[entry]['Cross-reference (EMBL)']):
                nucleicFilePath = '{}/{}.embl'.format(tmpDirectoryProcess, getNucleicFialeName(nucleicAccession))
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
                    outputContent.append([
                        entry,
                        crossReference[entry]['Cross-reference (embl)'][index],
                        'protein_id',
                        nucleicAccession,
                        'embl',
                        nucleicFilePath
                    ])
        with open(outputName, 'w') as file:
            file.write('{}\n'.format('\t'.join(common.global_dict['inputIIheaders'])))
            for line in ['\t'.join(values) for values in outputContent]:
                file.write('{}\n'.format(line))
        logger.info('{} generated.'.format(outputName))
        logger.info('{} finished!'.format(boxName))
    else:
        logger.error('Input header unrecognized.')
        exit(1)
