#!/usr/bin/env python3

import re
import sys
import json
import urllib3
import certifi
import xml.etree.ElementTree as ET
import subprocess
import os
import shutil
from Bio import SeqIO

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

def check_inputFile(lines):
    ''' controles the consistency of the input file provided by the user
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*- errors or warnings must be implemented in problematic     -*-*-*-
    -*-*-*- cases                                                     -*-*-*-
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    '''
    d_check = {}
    for aline in lines[1:]:
        # print(line[lines.index(aline)])
        splitted = aline.strip().split("\t")
        if splitted[0] in d_check:
            if splitted[2] != d_check[splitted[0]]:
                print("A protein accession cannot appear for different nucleotide accessions")
                print(lines[lines.index(aline)])
        else:
            d_check[splitted[0]] = splitted[2]
        if splitted[2] in d_check:
            if splitted[4] != d_check[splitted[2]][0]:
                print("A nucleotide accession cannot appear in different nucleotide files")
                print(lines[lines.index(aline)])
            if splitted[5] != d_check[splitted[2]][1]:
                if d_check[splitted[2]][1] == "NA":
                    d_check[splitted[2]] = [splitted[4], splitted[5]]
                else:
                    print("A nucleotide accession cannot have various taxonID")
                    print(lines[lines.index(aline)])
        else:
            d_check[splitted[2]] = [splitted[4], splitted[5]]
        if splitted[4] in d_check:
            if splitted[1] != d_check[splitted[4]][0]:
                print("The protein accession field must be consistent all along a nucleotide file")
                print(lines[lines.index(aline)])
            if splitted[3] != d_check[splitted[4]][1]:
                print("The nucleotide file format must be consistent all along a nucleotide file")
                print(lines[lines.index(aline)])
        else:
            d_check[splitted[4]] = [splitted[1], splitted[3]]
    return 0

def create_d_input(lines):
    ''' lines are taken from the input provided by the user
    the function creates a dictionary recording the information about the
    targets (protAC), the associated field (protACfield), the nucleotide
    accession (nucAC, = contig identifier), the INSDC file format
    (nucFileFormat), the file name (nucFile) and the taxonID

    If the taxonID is mentioned only once while several lines are related to
    the same file, it is able to manage with the lack of information, else, the
    taxonID must be provided inside the INSDC file. The user information has the
    priority over INSDC file taxonID
    '''
    d_input = {}
    for aline in lines[1:]: # skips the first line, with headers
        aline = aline.strip().split("\t")
        d_input.setdefault(aline[4], {}).setdefault(aline[2], {}).setdefault("target_list", [])
        if aline[0] not in d_input[aline[4]][aline[2]]["target_list"]:
            d_input[aline[4]][aline[2]]["target_list"].append(aline[0])
        d_input[aline[4]]["protACfield"] = aline[1]
        d_input[aline[4]]["nucFileFormat"] = aline[3]
        if "taxonID" in d_input[aline[4]][aline[2]]:
            if d_input[aline[4]][aline[2]]["taxonID"] == "NA":
                del d_input[aline[4]][aline[2]]["taxonID"]

        d_input[aline[4]][aline[2]].setdefault("taxonID", aline[5])
    return d_input

def check_and_get_input(input):
    ''' 
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    -*-*-*- faire le check en même temps que la création du dico -*-*-*-
    -*-*-*- revoir la liste des checks à faire -*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    '''
    with open(input, "r") as infile:
        lines = infile.readlines()
    lines = list(skip_duplicates(lines)) #remove duplicates

    check_inputFile(lines) #has to be rethought and integrate the check
    d_input = create_d_input(lines) #while creating input
    return d_input

def get_from_dbxref(aFeature, dbref):
    ''' retrieves the value from the dbxref list
    '''
    if dbref == "taxon":
        pattern = "taxon:"
    elif dbref == "MaGe":
        pattern = "MaGe:"
    elif dbref == "UniProt":
        pattern = "UniProt*"

    if aFeature.qualifiers.get('db_xref'):
        for aRef in aFeature.qualifiers.get('db_xref'):
            if re.match(pattern, aRef):
                return aRef.split(r':')[-1]
            else:
                return "NA"
    return "NA"

def get_uniq_value(aFeature, ref):
    ''' makes sure that when there is a value, it is not a sequence of values
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*- in that case a value is a sequence of values, it has to   -*-*-*-
    -*-*-*- send an error, or a warning (need to know which index to  -*-*-*-
    -*-*-*- select)                                                   -*-*-*-
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    '''
    result = aFeature.qualifiers.get(ref)
    if not result:
        try:
            result = get_from_dbxref(aFeature, ref)
            return result
        except:
            return "NA"
    if len(result) == 1:
        return result[0]
    else:
        print("warning !!!", result)
        return 1

def get_required_value(func, aFeature, *args):
    ''' function called if the value is mandatory (protein ID (ident),
    sequence)
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*- an error has to be implemented so that any mandatory value -*-*-*-
    -*-*-*- couldn't be None or "NA"                                   -*-*-*-
    -*-*-*- maybe do a try/except evaluation
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    '''
    if not func(aFeature, *args) or func(aFeature, *args) == "NA":
        print("there is a rho rho problem")
    else:
        return func(aFeature, *args)

def search_taxonID(aFeature, given_taxon):
    ''' makes sure that the taxon ID is provided by the user or through the
    INSDC file
    '''
    if given_taxon != "NA":
        #print("taxon_id provided by the user")
        taxon_id = given_taxon
    else:
        #print("taxon_id must be in the file")
        taxon_id = get_required_value(get_uniq_value, aFeature, "taxon")
    return taxon_id

def get_contig_info(aFeature, contig_content, given_taxon):
    ''' then get the relied information to a contig as organism, strain, size
    and taxon ID
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*- an error has to be implemented so that any cds couldn't be -*-*-*-
    -*-*-*- relied to a taxon ID                                       -*-*-*-
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    '''
    taxon_id = search_taxonID(aFeature, given_taxon)
    if taxon_id != "NA":
        contig_content.update({
            "organism": get_uniq_value(aFeature, "organism"),
            "strain": get_uniq_value(aFeature, "strain"),
            "taxon_id": taxon_id,
            "size": [aFeature.location.start.real+1,
                     aFeature.location.end.real
                    ],
            "cds_to_keep": [],
            "window": []
            })
    else:
        print("We have a problem")
    return contig_content

def is_pseudogene(aFeature):
    ''' tests if a CDS is a pseudogene
    '''
    if "pseudo" in aFeature.qualifiers:
        return True
    else:
        return False

def get_pseudo_id(params):
    ''' creates an identifier for pseudogenes
    '''
    params["INC_PSEUDO_REF"] += 1
    INC_PSEUDO_REF = params["INC_PSEUDO_REF"]
    return ":".join(["PSEUDO", str(INC_PSEUDO_REF)])

def det_frame(strand, startCDS, endGenome):
    ''' determines the frame of the CDS
    '''
    if strand == "1":
        window = ((int(startCDS)-1)%3)+1
        return "".join(["+", str(window)])
    elif strand == "-1":
        window = ((int(endGenome)-int(startCDS))%3)+1
        return "".join(["-", str(window)])

def get_pseudo_info(aFeature, cds_info, contig_content, params):
    ''' adds to the window's list information on the pseudogene
    '''
    params["INC_CDS_REF"] += 1
    INC_CDS_REF = params["INC_CDS_REF"]

    pseudo_id = get_pseudo_id(params)
    start, stop = aFeature.location.start.real+1, aFeature.location.end.real
    cds_info[INC_CDS_REF] = {
        "protein_id": pseudo_id,
        "position": [start, stop],
        "frame": det_frame(str(aFeature.location.strand), start, contig_content.get("size")[1]),
        "product": get_uniq_value(aFeature, "product"),
        "sequence": get_uniq_value(aFeature, "translation")
        }
    contig_content["window"].append((INC_CDS_REF,pseudo_id))
    return cds_info, contig_content, params

def get_prot_info(aFeature, cds_info, contig_content, proteinField, params):
    ''' adds to the window's list information on the protein
    '''
    INC_CONTIG_REF = params["INC_CONTIG_REF"]
    params["INC_CDS_REF"] += 1
    INC_CDS_REF = params["INC_CDS_REF"]

    ident = get_required_value(get_uniq_value, aFeature, proteinField) 
    start, stop = aFeature.location.start.real+1, aFeature.location.end.real
    ec_nums = (aFeature.qualifiers.get("EC_number") if aFeature.qualifiers.get("EC_number") else "NA")
    cds_info[INC_CDS_REF] = {
        "protein_id": ident,
        "position": [start, stop],
        "frame": det_frame(str(aFeature.location.strand), start, contig_content.get("size")[1]),
        "product": get_uniq_value(aFeature, "product"),
        "sequence": get_required_value(get_uniq_value, aFeature, "translation"),
        "ec_number": ec_nums,
        "uniprot": get_from_dbxref(aFeature, "UniProt"),
        "target": False,
        "context": None,
        "contig": INC_CONTIG_REF
        }
    contig_content["window"].append((INC_CDS_REF, ident))
    return cds_info, contig_content, params

def is_target(cds, target_list, cds_info, contig_content, params):
    ''' tests if the tested CDS is referenced as a target
    If yes, copy the window's information to cds_keeper :
    information on target and genomic context are saved

    *** tested value used to be the CDS in the middle of the window ***
    '''
    cds_ref = cds[0]
    cds_id = cds[1]
    if cds_id in target_list:
        #print(contig_content["window"])
        contig_content["cds_to_keep"].extend(contig_content["window"])
        #print("> ", contig_content["window"])
        cds_info[cds_ref].update({"target": True ,
                                 "context": contig_content["window"].copy()
                                 })
        #print(">> ", cds_info[cds_ref]["context"])
        params["INC_TARGET_LOADED"] += 1
    return cds_info, contig_content, params

def is_kept(cds, cds_info, contig_content):
    ''' verifies if the tested CDS is in cds_keeper list
    If not, all information on it are removed

    *** tested value is the first value of the window ***
    '''
    cds_ref = cds[0]
    cds_id = cds[1]
    if contig_content["cds_to_keep"]:
        if cds_id not in [ref[1] for ref in contig_content["cds_to_keep"]]:
            del cds_info[cds_ref]
    else:
        del cds_info[cds_ref]
    del contig_content["window"][0]
    return cds_info, contig_content

def sliding_window(cds_info, contig_content, target_list, beginContig, params):
    ''' makes the window going ahead through the CDSs in a contig
    '''
    MAX_GC = params["MAX_GC"]
    HALF_SIZE_GC = int(MAX_GC/2)
    window = contig_content["window"]
    window_size = len(window)
    if window_size == MAX_GC:
        beginContig = False
        cds_info, contig_content, params = is_target(window[HALF_SIZE_GC], target_list, cds_info, contig_content, params)
        cds_info, contig_content = is_kept(window[0], cds_info, contig_content)
    # COM: window's size didn't reach the MAX_GC yet
    elif window_size > HALF_SIZE_GC and beginContig:
        cds_info, contig_content, params = is_target(window[window_size-(HALF_SIZE_GC+1)], target_list, cds_info, contig_content, params)
    # COM: the window is overrunning the contig
    elif window_size > HALF_SIZE_GC and not beginContig:
        cds_info, contig_content, params = is_target(window[HALF_SIZE_GC], target_list, cds_info, contig_content, params)
        cds_info, contig_content = is_kept(window[0], cds_info, contig_content)
    # COM: there is no more target at the end of the contig
    elif window_size <= HALF_SIZE_GC and not beginContig:
        cds_info, contig_content = is_kept(window[0], cds_info, contig_content)
    return cds_info, contig_content, beginContig, params

def parse_insdc(afile, d_infile, cds_info, contig_info, params):
    STOP_INC_CONTIG = params["INC_CONTIG_REF"] + len(d_infile)-2

    typeParsing = d_infile["nucFileFormat"]
    fieldProteinID = d_infile["protACfield"]
    
    with open(afile, "r") as insdcFile:
        seqRecordParsed = SeqIO.parse(insdcFile, typeParsing)
        for seqRecord in seqRecordParsed:
            if params["INC_CONTIG_REF"] >= STOP_INC_CONTIG:
                break
            contig_name = seqRecord.id.split(r'.')[0]
            if contig_name in d_infile:
                params["INC_CONTIG_REF"] += 1
                INC_CONTIG_REF = params["INC_CONTIG_REF"]
                params["INC_TARGET_LOADED"] = 0
                beginContig = True
                contig_info[INC_CONTIG_REF] = {"contig": contig_name}
                for aFeature in seqRecord.features:
                    if aFeature.type == "source":
                        contig_info[INC_CONTIG_REF] = get_contig_info(
                            aFeature,
                            contig_info[INC_CONTIG_REF],
                            d_infile[contig_name]["taxonID"]
                            )
                    elif aFeature.type == "CDS":
                        newCdsAdded = False
                        if params["INC_TARGET_LOADED"] >= len(d_infile[contig_name]["target_list"]):
                            break
                        if is_pseudogene(aFeature):
                            if params["PSEUDOGENE"]:
                                cds_info, contig_info[INC_CONTIG_REF], params = get_pseudo_info(aFeature, cds_info, contig_info[INC_CONTIG_REF], params)
                                newCdsAdded = True
                            #else: # ???
                                #break # ???
                        else:
                            cds_info, contig_info[INC_CONTIG_REF], params = get_prot_info(aFeature, cds_info, contig_info[INC_CONTIG_REF], fieldProteinID, params)
                            # print(params, [contig_info[ref]["window"]
                            #                for ref in contig_info])                                           
                            newCdsAdded = True
                        if newCdsAdded:
                            cds_info, contig_info[INC_CONTIG_REF], beginContig, params = sliding_window(
                                cds_info,
                                contig_info[INC_CONTIG_REF],
                                d_infile[contig_name]["target_list"],
                                beginContig,
                                params
                                )
                for i in range(params["MAX_GC"]-1):# don't work with "_" instead of "i"
                    cds_info, contig_info[INC_CONTIG_REF], beginContig, params = sliding_window(
                        cds_info,
                        contig_info[INC_CONTIG_REF],
                        d_infile[contig_name]["target_list"],
                        beginContig,
                        params
                        )
                contig_info[INC_CONTIG_REF]["cds_to_keep"] = list(skip_duplicates(contig_info[INC_CONTIG_REF]["cds_to_keep"]))

    return cds_info, contig_info, params


def parse_INSDC_files(d_input, cds_info, contig_info, params):
    ''' 
    '''
    for afile in d_input:
        cds_info, contig_info, params = parse_insdc(afile, d_input[afile], cds_info, contig_info, params)
        #print(params["INC_CONTIG_REF"], len(d_input[afile])-2)
    return cds_info, contig_info, params

def concat_by_dot(prefix, suffix):
    ''' does the concatenation by a dot
    '''
    return ".".join([prefix, suffix])

def write_multiFasta(cds_info, prefix, suffix):
    ''' writes a multiFasta file from the dictionary obtained by the INSDC file
    parser
    '''
    fileToWrite = concat_by_dot(prefix, suffix)
    with open(fileToWrite, "w") as fastaFile:
        for cds in cds_info:
            count = 0
            fastaFile.write("".join([">", cds_info[cds]["protein_id"], "\n"]))
            for ac_amine in cds_info[cds]["sequence"]:
                count += 1
                fastaFile.write(ac_amine)
                if count >= 80:
                    fastaFile.write("\n")
                    count = 0
            fastaFile.write("\n")
    return 0

def write_json(dictionary, prefix, suffix):
    ''' writes all dictionaries of INSDC file parsed in a json file format
    '''
    filename = concat_by_dot(prefix, suffix)
    with open(filename, "w") as jsonFile:
        json.dump(dictionary, jsonFile, indent=2)
    return 0

def get_lineage(xml):
    ''' extracts the taxonomic lineage from the provided xml file format
    '''
    lineage_full = {}
    root = ET.fromstring(xml)
    scientificName = root.find('taxon').get('scientificName')
    rank = False
    if root.find('taxon').get('rank'):
        rank = root.find('taxon').get('rank')
    taxId = root.find('taxon').get('taxId')
    lineage_full[scientificName] = {'rank' : rank, 'taxId' : taxId}

    lineage = root.find('taxon').find("lineage")
    taxons = lineage.findall("taxon")
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
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    xml = http.request('GET', 'https://www.ebi.ac.uk/ena/data/view/Taxon:' + str(taxonID) + '&display=xml')
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
    if not desired_ranks == {}:
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
    for ataxon in taxonIDs:
        res_taxo = get_taxo_from_web(ataxon)
        desired_taxo = get_desired_lineage(res_taxo)
        new_lineage = store_into_dict(ataxon, desired_taxo)
        all_lineages.update(new_lineage) # no risk of overwriting
    return all_lineages

def mmseqs_createdb(prefix):
    ''' creates a database using the mmseqs software
    '''
    suffix = "faa"
    multiFasta = concat_by_dot(prefix, suffix)
    with open('mmseqs_createdb.log', 'w') as file:
        db_creation = subprocess.run(["mmseqs", "createdb", multiFasta, "ALL.DB"], stdout=file, stderr=file, check=True)
        print('createdb - exit code: {}'.format(db_creation.returncode))
    return 0

def mmseqs_clustering(prefix, clust_mode, cov, ident, cov_mode, cascaded):
    ''' does the clustering using the mmseqs software
    '''
    suffixDB = "DB"
    suffixCluster = "cluster"
    if os.path.isdir("MMseqsTMP"):
        shutil.rmtree("MMseqsTMP")
    os.mkdir("MMseqsTMP")
    dataBase = concat_by_dot(prefix, suffixDB)
    outputCluster = concat_by_dot(prefix, suffixCluster)
    with open('mmseqs_clustering.log', 'w') as file:
        clust_creation = subprocess.run(["mmseqs", "cluster", dataBase,
                        outputCluster, "MMseqsTMP",
                        "--min-seq-id", str(ident),
                        "--cov-mode", str(cov_mode),
                        "-c", str(cov),
                        "--cluster-mode", str(clust_mode)#,
                        #"--cascaded", str(cascaded)
                       ], stdout=file, stderr=file, check=True)
        print('clustering - exit code: {}'.format(clust_creation.returncode))
    return 0

def mmseqs_createTSV(prefix):
    ''' executes the mmseqs command line "mmseqs createtsv"
    '''
    suffixDB = "DB"
    suffixCluster = "cluster"
    suffixTSV = "tsv"
    inputDB = concat_by_dot(prefix, suffixDB)
    inputCluster = concat_by_dot(prefix, suffixCluster)
    outputTSV = concat_by_dot(prefix, suffixTSV)
    with open('mmseqs_createtsv.log', 'w') as file:
        tsv_creation = subprocess.run(["mmseqs", "createtsv", inputDB,
                        inputDB, inputCluster, outputTSV
                       ], stdout=file, stderr=file, check=True)
        print('createTSV - exit code: {}'.format(tsv_creation.returncode))
    return 0

def mmseqs_runner(params):
    ''' runs the mmseqs2 software on the multiFasta file "ALL.faa"
    '''
    mmseqs_createdb(params["prefix"])
    mmseqs_clustering(params["prefix"], params["cluster_mode"], params["coverage"],
                      params["min_id"], params["cov_mode"], params["cascaded"])
    mmseqs_createTSV(params["prefix"])
    return 0

def regroup_families(prefix):
    ''' creates a dictionary to store families obtained by MMseqs2
    '''
    suffix = "tsv"
    INC_FAMILY = 1
    family_in_progress = []
    families = {}
    with open(concat_by_dot(prefix, suffix), "r") as file:
        lines = file.readlines()
    for aline in lines:
        aline = aline.strip().split("\t")
        if aline[0] in family_in_progress:
            family_in_progress.append(aline[1])
        else:
            fam_name = "_".join(["family", str(INC_FAMILY)])
            if family_in_progress != []:
                families[fam_name] = family_in_progress
                INC_FAMILY += 1
            family_in_progress = [aline[1]]
    fam_name = "_".join(["family", str(INC_FAMILY)])
    families[fam_name] = family_in_progress
    return families

def run(input_file, prefix, params):
    ''' main script to run the second box of NetSyn2
    '''

    params = {
        "PSEUDOGENE": False, # Tells if pseudogenes are included in the analysis
        "MAX_GC": 5, # size of the window
        "INC_PSEUDO_REF": 0, # counter of pseusogenes
        "INC_CDS_REF": 0,
        "INC_CONTIG_REF": 0,
        "min_id": 0.3,
        "cov_mode": 1,
        "coverage": 0.8,
        "cluster_mode": 0,
        "cascaded": False
        }

    try:
        shutil.rmtree("MMseqsTMP/")
    #print(os.listdir())
        os.remove("ALL.cluster")
        os.remove("ALL.cluster.index")
        os.remove("ALL.tsv")
        os.remove("mmseqs_createtsv.log")
    except:
        pass


    params["prefix"] = prefix
    cds_info = {}
    contig_info = {}
    d_input = check_and_get_input(input_file)
    #tester la fonction map() de python pour appliquer une fonction sur une
    #liste
    #usage : map(myFun, myList)
    cds_info, contig_info, params = parse_INSDC_files(d_input, cds_info, contig_info, params)

    # for akey in cds_info:
    #     print(akey, cds_info[akey], sep='\t')
    #     print("***")
    # print('\n~~~~~~~~~~~~~~\n')
    # for akey in contig_info:
    #     print(akey, contig_info[akey], sep='\t')
    #     print("***")

    write_multiFasta(cds_info, prefix, "faa")
    write_json(cds_info, prefix, "GC.json")
    write_json(contig_info, prefix, "contig.json")
    
    taxonIDs = list(set([contig_info[contig]["taxon_id"] for contig in contig_info]))
    taxonomicLineage = get_taxonLineage(taxonIDs)
    write_json(taxonomicLineage, prefix, "lineage.json")

    mmseqs_runner(params)
    
    families = regroup_families(prefix)
    real_families = {key: value for (key, value) in families.items() if len(value) > 1}
    singletons = {key: value for (key, value) in families.items() if len(value) == 1}
    print(real_families)
    write_json(real_families, prefix, "families.json")
    write_json(singletons, prefix, "singletons.json")

if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2])
