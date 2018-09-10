import xml.etree.ElementTree as ET
import re
import os
import shutil
import sys
import json
import subprocess
import urllib3
import certifi
from Bio import SeqIO

# if FORMAT_INSDC == "ncbi_embl":
#     afile = "ncbi_files/acineto.embl"
#     d_input = {
#         "CU468230":["CAO99451.1", "CAO99632.1", "CAO99422.1"]
#         }
# elif FORMAT_INSDC == "mic_embl":
#     afile = "mic_files/ABSDF.1.embl"
#     d_input = {
#         "ABSDF":["2392187", "2392438", "2392146", "2392153"]
#         }
# elif FORMAT_INSDC == "ncbi_gbk":
#     afile = "ncbi_files/acineto.gbk"
#     d_input = {
#         "CU468230":["CAO99426.1", "CAO99587.1", "CAO99702.1"]
#         }
# elif FORMAT_INSDC == "mic_gbk":
#     afile = "mic_files/ABSDF.1.gbk"
#     d_input = {
#         "ABSDF":["2392146", "2392228", "2392230", "2395821"]
#         }
# elif FORMAT_INSDC == "ncbi_gbff":
#     afile = "ncbi_files/acineto.gbff"
#     d_input = {
#         "CU468230":["CAO99451.1", "CAO99742.1"]
#         }
# elif FORMAT_INSDC == "mic_gbff":
#     afile = "mic_files/NZ_ADHA.1-Contigs.gbff"
#     d_input = {
#         "NZ_ADHAC00001":["5537943"],
#         "NZ_ADHAC00002":["5537948"],
#         "NZ_ADHAC00005":["5538000", "5538003"]
#         }
# elif FORMAT_INSDC == "gbff":
#     contig_list = ["CU468230"]
#     target_list = 


PSEUDOGENE = False # Tells if pseudogenes are included in the analysis
MAX_GC = 5 # size of the window
INC_PSEUDO = 0 # counter of pseusogenes
IN_CONTIG = False
IN_PSEUDO = False
parserINSDCfile = {}

def parse_info(d_file):
    ''' returns information about the file to parse it correctly
    '''
    global TYPE_PARSING, FIELD_PROTEIN_ID
    TYPE_PARSING = d_file["nucFileFormat"]
    FIELD_PROTEIN_ID = d_file["protACfield"]
    return TYPE_PARSING, FIELD_PROTEIN_ID

def det_frame(strand, startCDS, endGenome):
    ''' determines the frame of the CDS
    '''
    if strand == "1":
        window = ((int(startCDS)-1)%3)+1
        return "".join(["+", str(window)])
    elif strand == "-1":
        window = ((int(endGenome)-int(startCDS))%3)+1
        return "".join(["-", str(window)])

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

def search_for_taxonID(given_taxon, aFeature):
    ''' makes sure that the taxon ID is provided by the user or through the
    INSDC file
    '''
    if given_taxon != "NA":
        print("provided")
        taxon_id = given_taxon
    else:
        print("must be in the file")
        taxon_id = get_from_dbxref(aFeature, "taxon")
    return taxon_id

def get_contig_info(aFeature, parser_contig_content, given_taxon):
    ''' then get the relied information to a contig as organism, strain, size
    and taxon ID
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*- an error has to be implemented so that any cds couldn't be -*-*-*-
    -*-*-*- relied to a taxon ID                                       -*-*-*-
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    '''
    taxon_id = search_for_taxonID(given_taxon, aFeature)
    if taxon_id != "NA":
        parser_contig_content = {
            "organism": get_uniq_value(aFeature, "organism"),
            "strain": get_uniq_value(aFeature, "strain"),
            "taxon_id": taxon_id,
            "size": [aFeature.location.start.real+1,
                     aFeature.location.end.real
                    ]
            }
    else:
        print("We have a problem")
    return parser_contig_content

def get_required_value(func, aFeature, *args):
    ''' function called if the value is mandatory (protein ID (ident),
    sequence)
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    -*-*-*- an error has to be implemented so that any mandatory value -*-*-*-
    -*-*-*- couldn't be None or "NA"                                   -*-*-*-
    -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    '''
    if not func(aFeature, *args) or func(aFeature, *args) == "NA":
        print("there is a rho rho problem")
    else:
        return func(aFeature, *args)

def is_pseudogene(aFeature):
    ''' tests if a CDS is a pseudogene
    '''
    if "pseudo" in aFeature.qualifiers:
        return True
    else:
        return False

def get_pseudo_id():
    ''' creates an identifier for pseudogenes
    '''
    global INC_PSEUDO
    INC_PSEUDO += 1
    return ":".join(["PSEUDO", str(INC_PSEUDO)])

def get_pseudo_info(aFeature, parser_contig_content):
    ''' adds to the window's list information on the pseudogene
    '''
    pseudo_id = get_pseudo_id()
    start, stop = aFeature.location.start.real+1, aFeature.location.end.real
    parser_contig_content[pseudo_id] = {
        "position": [start, stop],
        "frame": det_frame(str(aFeature.location.strand), start, parser_contig_content.get("size")[1]),
        "product": get_uniq_value(aFeature, "product"),
        "sequence": get_uniq_value(aFeature, "translation")
        }
    parser_contig_content["window"].append(pseudo_id)
    return parser_contig_content

def get_prot_info(aFeature, parser_contig_content):
    ''' adds to the window's list information on the protein
    '''
    ident = get_required_value(get_uniq_value, aFeature, FIELD_PROTEIN_ID)
    parser_contig_content.setdefault(ident, {})
    start, stop = aFeature.location.start.real+1, aFeature.location.end.real
    ec_nums = (aFeature.qualifiers.get("EC_number") if aFeature.qualifiers.get("EC_number") else "NA")
    parser_contig_content[ident] = {
        "position": [start, stop],
        "frame": det_frame(str(aFeature.location.strand), start, parser_contig_content.get("size")[1]),
        "product": get_uniq_value(aFeature, "product"),
        "sequence": get_required_value(get_uniq_value, aFeature, "translation"),
        "ec_number": ec_nums,
        "uniprot": get_from_dbxref(aFeature, "UniProt")
        }
    parser_contig_content["window"].append(ident)
    return parser_contig_content

def is_target(tested, target_list, window, cds_keeper):
    ''' tests if the tested CDS is referenced as a target
    If yes, copy the window's information to cds_keeper :
    information on target and genomic context are saved

    *** tested value used to be the CDS in the middle of the window ***
    '''
    if tested in target_list:
        cds_keeper.extend(window)
        global INC_TARGETS_LOADED
        INC_TARGETS_LOADED += 1
    return cds_keeper

def is_kept(tested, cds_keeper, parser_contig_content, window):
    ''' verifies if the tested CDS is in cds_keeper list
    If not, all information on it are removed

    *** tested value is the first value of the window ***
    '''
    if cds_keeper:
        if tested not in cds_keeper:
            del parser_contig_content[tested]
            #return parser_contig_content
    else:
        del parser_contig_content[tested]
    del window[0]
    return parser_contig_content, window

def sliding_window(parser_contig_content, BEGIN_CONTIG, MAX_GC, target_list):
    ''' makes the window going ahead through the CDSs in a contig
    '''
    window = parser_contig_content.get('window')
    cds_keeper = parser_contig_content.get('cds_to_keep')
    HALF_SIZE_GC = int(MAX_GC/2)
    window_size = len(window)
    if window_size == MAX_GC:
        BEGIN_CONTIG = False
        cds_keeper = is_target(window[HALF_SIZE_GC], target_list, window, cds_keeper)
        parser_contig_content, window = is_kept(window[0], cds_keeper, parser_contig_content, window)
    # COM: window's size didn't reach the MAX_GC yet
    elif window_size > HALF_SIZE_GC and BEGIN_CONTIG:
        cds_keeper = is_target(window[window_size-(HALF_SIZE_GC+1)], target_list, window, cds_keeper)
    # COM: the window is overrunning the contig
    elif window_size > HALF_SIZE_GC and not BEGIN_CONTIG:
        cds_keeper = is_target(window[HALF_SIZE_GC], target_list, window, cds_keeper)
        parser_contig_content, window = is_kept(window[0], cds_keeper, parser_contig_content, window)
    # COM: there is no more target at the end of the contig
    elif window_size <= HALF_SIZE_GC and not BEGIN_CONTIG:
        parser_contig_content, window = is_kept(window[0], cds_keeper, parser_contig_content, window)
    return parser_contig_content, BEGIN_CONTIG

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

def parseAfile(afile, d_afile):
    ''' does the parsing of an INSDC file
    Calls the functions in order to get information about contig, protein
    and pseudogene (if required by the user), and the sliding_window function
    All kept information are stored in a dictionary called parserINSDCfile
    '''
    TYPE_PARSING, FIELD_PROTEIN_ID = parse_info(d_afile)
    with open(afile) as input_file:
        print(afile)
        seqRecordParsed = SeqIO.parse(input_file, TYPE_PARSING)
        INC_CONTIG_LOADED = 0
        for seqRecord in seqRecordParsed:
            if INC_CONTIG_LOADED >= len(d_afile)-2:
                break
            contig_name = seqRecord.id.split(r'.')[0]
            if contig_name in d_afile.keys():
                INC_CONTIG_LOADED += 1
                global INC_TARGETS_LOADED
                INC_TARGETS_LOADED = 0 # counter of loaded targets
                BEGIN_CONTIG = True
                parserINSDCfile.setdefault(contig_name, {})
                for aFeature in seqRecord.features:
                    if aFeature.type == "source": #INFO_CONTIG:
                        print(d_afile[contig_name]["taxonID"])
                        parserINSDCfile[contig_name] = get_contig_info(aFeature, parserINSDCfile[contig_name], d_afile[contig_name]["taxonID"])
                        parserINSDCfile[contig_name].setdefault("cds_to_keep", [])
                        parserINSDCfile[contig_name].setdefault("window", [])
                    elif aFeature.type == "CDS":
                        NEW_CDS_ADDED = False
                        if INC_TARGETS_LOADED >= len(d_afile[contig_name]["target_list"]):
                            break
                        if is_pseudogene(aFeature):
                            if PSEUDOGENE: # defined by the user : takes into account pseudogenes
                                parserINSDCfile[contig_name] = get_pseudo_info(aFeature, parserINSDCfile[contig_name])
                                NEW_CDS_ADDED = True
                        else:
                            parserINSDCfile[contig_name] = get_prot_info(aFeature, parserINSDCfile[contig_name])
                            NEW_CDS_ADDED = True
                        if NEW_CDS_ADDED:
                            parserINSDCfile[contig_name], BEGIN_CONTIG = sliding_window(
                                parserINSDCfile[contig_name],
                                BEGIN_CONTIG,
                                MAX_GC,
                                d_afile[contig_name]["target_list"])
                for i in range(MAX_GC-1):
                    parserINSDCfile[contig_name], BEGIN_CONTIG = sliding_window(
                        parserINSDCfile[contig_name],
                        BEGIN_CONTIG,
                        MAX_GC,
                        d_afile[contig_name]["target_list"])
                parserINSDCfile[contig_name]["cds_to_keep"] = list(skip_duplicates(parserINSDCfile[contig_name]["cds_to_keep"]))
    return parserINSDCfile

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
        d_input.setdefault(aline[4], {}).setdefault(aline[2], {}).setdefault("target_list", []).append(aline[0])
        d_input[aline[4]]["protACfield"] = aline[1]
        d_input[aline[4]]["nucFileFormat"] = aline[3]

        if "taxonID" in d_input[aline[4]][aline[2]]:
            if d_input[aline[4]][aline[2]]["taxonID"] == "NA":
                del d_input[aline[4]][aline[2]]["taxonID"]

        d_input[aline[4]][aline[2]].setdefault("taxonID", aline[5])
    return d_input

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
                print("A nucleotide accession caonnot appear in different nucleotide files")
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

def write_multiFasta(aFileParsed):
    ''' writes a multiFasta file from the dictionary obtained by the INSDC file
    parser
    '''
    with open("ALL.faa", "w") as fastaFile:
        for acont in aFileParsed:
            for acds in aFileParsed[acont]["cds_to_keep"]:
                count = 0
                fastaFile.write("".join([">", acds, "\n"]))
                for aaa in aFileParsed[acont][acds]['sequence']:
                    count += 1
                    fastaFile.write(aaa)
                    if count >= 80:
                        fastaFile.write("\n")
                        count = 0
                fastaFile.write("\n")
    return 0

def write_json(aFileParsed, filename):
    ''' writes all dictionaries of INSDC file parsed in a json file format
    '''
    with open(filename, "w") as jsonFile:
        json.dump(aFileParsed, jsonFile, indent=4)
    return 0

def get_taxonIDs(aFileParsed):
    ''' get a non redundant list of taxon IDs represented in the analysis
    '''
    taxon_list = []
    for acont in aFileParsed:
        taxon_list.append(aFileParsed[acont]["taxon_id"])
    return list(set(taxon_list))

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

# def turn_into_dico(ataxon, desired_taxo, all_lineage):
#     d_alineage = {}
#     desired_taxo.sort()
#     print(desired_taxo, "\n")
#     IS_SON = False
#     for alist_info in desired_taxo:
#         if alist_info[3] != "NA":
#             if alist_info[0] == 1:
#                 father = alist_info[3]
#             elif not father:
#                 father = "NA"
#             if IS_SON:
#                 d_alineage[father]["descendants"].append(alist_info[3])
#                 IS_SON = False
#             d_alineage[alist_info[3]] = {
#                 "scientificName": alist_info[2],
#                 "level": alist_info[0],
#                 "rank": alist_info[1],
#                 "father": father,
#                 "descendants": []
#                 }
#             if desired_taxo.index(alist_info) == len(desired_taxo)-1:
#                 if ataxon != alist_info[3]:
#                     d_alineage[alist_info[3]]["descendants"].append(ataxon)
#             father = alist_info[3]
#         else:
#             IS_SON = True
#     return d_alineage
#
# def check_new_lineage(all_lineages, new_lineage):
#     for new_key in new_lineage:
#         if new_key in all_lineages:
#             if new_lineage[new_key] != all_lineages[new_key]:
#                 for key, value in new_lineage[new_key].items():
#                     if value != all_lineages[new_key][key]:
#                         if key != "descendants":
#                             print("There are divergent lineage taxonomy or gap that might be filled")
#     return new_lineage

def store_into_dico(ataxon, desired_taxo):
    ''' puts information on taxonomic lineage into a dictionary
    '''
    d_newLineage = {}
    d_newLineage[ataxon] = {}
    for alevel in desired_taxo:
        d_newLineage[ataxon][alevel[1]] = [level_info for level_info in alevel
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
        new_lineage = store_into_dico(ataxon, desired_taxo)
        all_lineages.update(new_lineage) # no risk of overwriting
    return all_lineages

def concat_by_dot(prefix, suffix):
    ''' does the concatenation by a dot
    '''
    return ".".join([prefix, suffix])

def mmseqs_createdb(prefix):
    ''' creates a database using the mmseqs software
    '''
    suffix = "faa"
    multiFasta = concat_by_dot(prefix, suffix)
    with open('mmseqs_createdb.log', 'w') as file:
        db_creation = subprocess.run(["mmseqs", "createdb", multiFasta, "ALL.DB"], stdout=file, stderr=file)
        print('exit code: {}'.format(db_creation.returncode))
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
        subprocess.run(["mmseqs", "cluster", dataBase,
                        outputCluster, "MMseqsTMP",
                        "--min-seq-id", str(ident),
                        "--cov-mode", str(cov_mode),
                        "-c", str(cov),
                        "--cluster-mode", str(clust_mode)#,
                        #"--cascaded", str(cascaded)
                       ], stdout=file, stderr=file)
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
        subprocess.run(["mmseqs", "createtsv", inputDB,
                        inputDB, inputCluster, outputTSV
                       ], stdout=file, stderr=file)
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

def main(input_file):
    ''' main script to run the second box of NetSyn2
    '''
    with open(input_file) as infile:
        lines = infile.readlines()
    lines = list(skip_duplicates(lines)) #remove duplicates

    check_inputFile(lines)
    d_input = create_d_input(lines)
    for afile in d_input:
        aFileParsed = parseAfile(afile, d_input[afile])

    write_multiFasta(aFileParsed)
    write_json(aFileParsed, "GC.json")

    taxonIDS_list = get_taxonIDs(aFileParsed)
    taxonomicLineage = get_taxonLineage(taxonIDS_list)
    write_json(taxonomicLineage, "TaxonomyLineage.json") # write TaxonomyLineage.json file

    params = {
        "prefix": "ALL",
        "min_id": 0.3,
        "cov_mode": 1,
        "coverage": 0.8,
        "cluster_mode": 0,
        "cascaded": False
        }
    mmseqs_runner(params)

    families = regroup_families(params["prefix"])
    write_json(families, "Families.json")

    print("END")

if __name__ == "__main__":
    main(sys.argv[1])
