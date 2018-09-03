##########
# Import #
##########
import os
import shutil
import logging
import urllib3
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


def getENAidCorrespondingToUniProtid(uniprotAccessions, batchesSize, PoolManager):
    '''
    Allows the correspondence between a UniProt accession and nuclotide accessions.
    Batch splitting to lighten the request.
    '''
    logger = logging.getLogger('{}.{}'.format(getENAidCorrespondingToUniProtid.__module__, getENAidCorrespondingToUniProtid.__name__))
    while uniprotAccessions:
        accessions = '+OR+id:'.join(uniprotAccessions[:batchesSize])
        try:
            res = PoolManager.request('GET' ,
                'https://www.uniprot.org/uniprot/?query=id:{}&columns=id,database(EMBL),database(EMBL_CDS)&format=tab'.format(accessions))
        except urllib3.exceptions:
            logger.error('OUPS')
            exit(1)
        print(res.data.decode('utf-8'))
        del uniprotAccessions[:batchesSize]

# def getEMBLfromENA(nucleotideAccession, PoolManager, logger):
    # res = PoolManager.request('GET' , 'https://www.ebi.ac.uk/ena/data/view/{}&display=text&set=true'.format(nucleotideAccession))

def run(InputName, tmpDirectory, maxGCsize, logger):
    '''
    Get INSDC files porocessing.
    '''
    boxName = 'GetINSDCFiles'
    print(dir(run))
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    logger.info('{} running...'.format(boxName))
    tmpDirectory = '{}/{}ProcessFiles'.format(tmpDirectory, boxName)
    if os.path.isdir(tmpDirectory):
        shutil.rmtree(tmpDirectory)
    os.mkdir(tmpDirectory)
    header, accessions = parseInputI(InputName)
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    http = urllib3.PoolManager()
    if header == 'UniProtAC':
        getENAidCorrespondingToUniProtid(accessions, 10, http)
    else:
        logger.error('Input header unrecognized.')
        exit(1)