#!/usr/bin/env python3

##########
# Import #
##########
from netsyn import common

import os
import logging
import gzip
import re
import time
import urllib3
import argparse
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry

# constants for uniprot query modification date 12/07/2022 
# constant for uniprot query
POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

#############
# Functions #
#############

def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    
    request.raise_for_status()
    session.close()
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        request.raise_for_status()
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
#                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(request["jobStatus"])
        else:
            session.close()
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    request.raise_for_status()
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
#    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    request.raise_for_status()
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results

def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/stream/")
    request = session.get(url)
    request.raise_for_status()
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)

def resultsFormat(res, dico):
    '''
    Formating ID match and update dico.
    '''
    logger = logging.getLogger('{}.{}'.format(
        resultsFormat.__module__, resultsFormat.__name__))
    sepColumn = '\t'
    sepField = ';'
    isHeader = True
    for bLine in res.data.splitlines():
        if isHeader:
            headers = bLine.decode('utf-8').split(sepColumn)
            if not headers[0] == 'Entry':
                logger.critical(headers)
                logger.critical('"id" column must be in first position')
                exit(2)
            isHeader = False
        else:
            for index, column in enumerate(bLine.decode('utf-8').split(sepColumn)):
                if index == 0:
                    entry = column
                    dico[entry] = {}
                elif not column:
                    logger.warning(
                        '{}: This entry is obsolete or not found. Please check in www.UniProt.org.'.format(entry))
                    del dico[entry]
                    break
                else:
                    dico[entry][headers[index]] = column.strip(
                        sepField).split(sepField)
    return dico


def getENAidMatchingToUniProtid(uniprotAccessions, batchesSize, PoolManager):
    '''
    Allows the correspondence between a UniProt accession and nuclotide accessions.9
    Batch splitting to lighten the request.
    '''
    logger = logging.getLogger('{}.{}'.format(
        getENAidMatchingToUniProtid.__module__, getENAidMatchingToUniProtid.__name__))
    crossReference = {}
    nbTotalEntries = len(uniprotAccessions)
    nbEntriesProcessed = 0
    logger.info(
        'Beginning of the correspondence between the UniProt and ENA identifiers...')
    crossReference = {}
    while uniprotAccessions:
        
#        accessions = '+OR+accession_id:'.join(uniprotAccessions[:batchesSize])
#        print("accesions {}".format(accessions))
#        res = common.httpRequest(
#            PoolManager, 'GET', 'https://www.uniprot.org/uniprot/?query=id:{}&columns=id,database(EMBL),database(EMBL_CDS)&format=tab'.format(accessions))
#            PoolManager, 'GET', 'https://rest.uniprot.org/uniprotkb/search?query=accession_id:{}&fields=accession,xref_embl&format=tsv'.format(accessions))
#        crossReference = resultsFormat(res, crossReference)
#        logger.info('uniprot {} crossref {}'.format(accessions, crossReference))
#        nbEntriesProcessed += len(uniprotAccessions[:batchesSize])

        accessions=uniprotAccessions[:batchesSize]
        job_id = submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="EMBL-GenBank-DDBJ_CDS", ids=accessions)
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)
            session.close()
            data=results['results']
            for row in data:
                if row['from'] not in crossReference:
                    crossReference[row['from']]={}
                    crossReference[row['from']]['Cross-reference (embl)']=[]
                    crossReference[row['from']]['Cross-reference (EMBL)']=[]
                crossReference[row['from']]['Cross-reference (embl)'].append(row['to'])

        job_id = submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="EMBL-GenBank-DDBJ", ids=accessions)
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)
            session.close()
            data=results['results']
            for row in data:
                 if row['from'] in crossReference:
                     crossReference[row['from']]['Cross-reference (EMBL)'].append(row['to'])

        del uniprotAccessions[:batchesSize]
        logger.info(
            'Correspondences computed: {}/{}'.format(nbEntriesProcessed, nbTotalEntries))
    return crossReference


def getNucleicFileName(nucleicAccession):
    '''
    Get nucleic file name from nucleic accession.
    '''
    nucleicAcc = re.match(
        r'^(?P<NucleicFileName>[A-Z]{4,6}[0-9]{2})[0-9]{6,}$', nucleicAccession)
    if nucleicAcc:
        return nucleicAcc.group("NucleicFileName")
    return nucleicAccession


def getEMBLfromENA(nucleicAccession, nucleicFilePath, PoolManager):
    '''
    Download EMBL file from ENA website.
    '''
    logger = logging.getLogger('{}.{}'.format(
        getEMBLfromENA.__module__, getEMBLfromENA.__name__))
    trial = True
    consTrialNB = 3
    trialNb = consTrialNB
    while trial:
        if trialNb < consTrialNB:
            time.sleep(2)
        trialNb -= 1
        res = common.httpRequest(
#            PoolManager, 'GET', 'https://www.ebi.ac.uk/ena/data/view/{}&display=text&set=true'.format(nucleicAccession))
#            PoolManager, 'GET', 'https://www.ebi.ac.uk/ena/browser/api/text/{}?lineLimit=0&annotationOnly=false&set=true'.format(nucleicAccession))
            PoolManager, 'GET', 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=genbank&retmode=text'.format(nucleicAccession)) 
        contentType = res.info()['Content-Type']
        print("contentType {}".format(contentType))
#        print("https://www.ebi.ac.uk/ena/browser/api/text/{}?lineLimit=0&annotationOnly=false&set=true".format(nucleicAccession))
        if contentType == 'text/plain;charset=UTF-8' and res.data.decode('utf-8') == 'Entry: {} display type is either not supported or entry is not found.\n'.format(nucleicAccession):
            logger.error(res.data.decode('utf-8'))
#            exit(1)
            break
        with open(nucleicFilePath, 'w') as file:
            if contentType == 'text/plain;charset=UTF-8':
                file.write(res.data.decode('utf-8'))
                trial = False
            elif contentType in ['application/x-gzip', 'application/octet-stream']:
                file.write(gzip.decompress(res.data).decode('utf-8'))
                trial = False
            elif contentType == 'text/plain':
                file.write(res.data.decode('utf-8'))
                trial = False
            else:
                if trialNb == 0 :
                    logger.critical(
                        'Unsupported content type ({}).'.format(contentType))
                    break
#                    trial = False
#                    exit(1)
    logger.info('{} downloaded.'.format(nucleicFilePath))


def run(InputName):
    '''
    Get INSDC files processing.
    '''
    # Constants
    boxName = common.global_dict['boxName']['GetINSDCFiles']
    dataDirectoryProcess = os.path.join(
        common.global_dict['dataDirectory'], boxName)
    outputName = common.global_dict['files'][boxName]['inputClusteringStep']
    # Logger
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    reportingMessages = []
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
        crossReference = getENAidMatchingToUniProtid(
            list(accessions), 200, http)
        withoutENAidNb = len(accessions)-len(crossReference)
        reportingMessages.append(
            'Targets without GENBANK correspondence number: {}/{}'.format(withoutENAidNb, len(accessions)))
        if withoutENAidNb:
            reportingMessages.append('Targets without GENBANK correspondence: {}'.format(
                set(accessions).difference(set(crossReference.keys()))
            ))
        outputContent = []
        logger.info('Beginning of GENBANK file downloading...')
        for entry in crossReference:
            maxAssemblyLength = 0
            for index, nucleicAccession in enumerate(crossReference[entry]['Cross-reference (EMBL)']):
                filesExtension = common.global_dict['filesExtension']

                # Previously, we were using the function `getNucleicFileName` as shown below:
                # nucleicFilePath = f'{os.path.join(dataDirectoryProcess, getNucleicFileName(nucleicAccession))}.{filesExtension}'
                # However, this caused issues because `getNucleicFileName` was removing too much from the accession,
                # resulting in incorrect file paths. We decided to use the accession directly instead.
                # This change can be reverted if needed in the future.

                nucleicFilePath = f'{os.path.join(dataDirectoryProcess, nucleicAccession)}.{filesExtension}'
                trial = True
                consTrialNb = 2
                trialNb = consTrialNb
                while trial:
                    if trialNb < consTrialNb:
                        os.remove(nucleicFilePath)
                        time.sleep(2)
                    trialNb -= 1
                    if not os.path.isfile(nucleicFilePath):
                        getEMBLfromENA(nucleicAccession, nucleicFilePath, http)
                    else:
                        logger.info(
                            '{} already existing.'.format(nucleicFilePath))
                    with open(nucleicFilePath, 'r') as file:
                        m = re.search(
#                            r'ID\s+{};.*; (?P<length>[0-9]+) BP\.\n'.format(nucleicAccession), file.read())
                            r'LOCUS\s+{}\s+ (?P<length>[0-9]+) bp'.format(nucleicAccession), file.read())  
                        if m:
                            assemblyLength = int(m.group('length'))
                            trial = False
                        else:
                            if not trialNb:
                                logger.error(
                                    '{}: improper modif file format.'.format(nucleicFilePath))
                                assemblyLength = 0
                                toAppend = []
                                break
#                                exit(1)
                if assemblyLength > maxAssemblyLength:
                    maxAssemblyLength = assemblyLength
                    for index2, nucleicAccession2 in enumerate(crossReference[entry]['Cross-reference (embl)']):
                        tmp_file = open(nucleicFilePath, 'r')
                        if nucleicAccession2 in tmp_file.read():
                            toAppend = [
                                entry,
                                crossReference[entry]['Cross-reference (embl)'][index2],
                                'protein_id',
                                nucleicAccession,
                                'genbank',
                                nucleicFilePath
                            ]
                            tmp_file.close()
                            break
                        tmp_file.close()
            if not toAppend:
                none = 0
            else:
                if not toAppend in outputContent:
                    outputContent.append(toAppend)
        withoutENAfilesNb = len(crossReference)-len(outputContent)
        reportingMessages.append(
            'Targets without EMBL file number: {}/{}'.format(withoutENAfilesNb, len(accessions)))
        if withoutENAfilesNb:
            reportingMessages.append('Targets without GENBANK file: {}'.format(
                set(crossReference.keys()).difference(
                    set([line[0] for line in outputContent]))
            ))
        reportingMessages.append('EMBL files number/organisms number downloded: {}'.format(
            len(set([line[5] for line in outputContent]))
        ))
        with open(outputName, 'w') as file:
            file.write('{}\n'.format(
                '\t'.join(common.global_dict['inputIIheaders'])))
            for line in ['\t'.join(values) for values in outputContent]:
                file.write('{}\n'.format(line))
        logger.info('{} generated.'.format(outputName))
        logger.info('{} completed!'.format(boxName))
        reportingMessages.append('Targets number conserved at end this step: {}/{}'.format(
            len(accessions)-withoutENAidNb-withoutENAfilesNb, len(accessions))
        )
        common.reportingFormat(logger, boxName, reportingMessages)
    else:
        logger.error('Input header unrecognized.')
#        exit(1)


def argumentsParser():
    '''
    Arguments parsing
    '''
    parser = argparse.ArgumentParser(description='version: {}'.format(common.global_dict['version']),
                                     usage='''GetINSDCFiles.py -u <UniProtAC.list> -o <OutputName>''',
                                     formatter_class=argparse.RawTextHelpFormatter)

    group1 = parser.add_argument_group('General settings')
    group1.add_argument('-u', '--UniProtACList', type=str,
                        required=True, help='Protein accession list.')
    group1.add_argument('-o', '--OutputName', type=str,
                        required=True, help='Name of corresponding file.')

    group2 = parser.add_argument_group('logger')
    group2.add_argument('--logLevel',
                        type=str,
                        nargs='?',
                        default='INFO',
                        help='log level',
                        choices=['ERROR', 'error', 'WARNING',
                                 'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                        required=False)
    group2.add_argument('--logFile',
                        type=str,
                        nargs='?',
                        help='log file (use the stderr by default)',
                        required=False)
    return parser.parse_args()


def main():
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
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault(
        'inputClusteringStep', '{}_correspondences.tsv'.format(args.OutputName))
    common.global_dict.setdefault('files', {}).setdefault(boxName, {}).setdefault(
        'report', '{}_{}_report.txt'.format(args.OutputName, boxName))

    #######
    # Run #
    #######
    run(args.UniProtACList)


if __name__ == '__main__':
    main()
