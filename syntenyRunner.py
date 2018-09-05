import sys
import json

GC_SIZE_USER = 3
HALF_SIZE_USER = int(GC_SIZE_USER/2)

def get_target_list(lines):
    target_list = []
    for aline in lines[1:]:
        aline = aline.strip().split("\t")
        target_list.append(aline[0])
    return list(set(target_list))

def find_index(target_in_contig, cds_list):
    indices = []
    for atarget in target_in_contig:
        indices.append(cds_list.index(atarget))
    return indices

def select_user_list(target_list, contig_dict):
    gc_user_list = []
    target_cds = {}
    target_in_contig = list(set(target_list) & set(contig_dict["cds_to_keep"]))
    indices = find_index(target_in_contig, contig_dict["cds_to_keep"])
    for index in indices:
        gc_chunk = contig_dict["cds_to_keep"][
            max(index-HALF_SIZE_USER, 0): min(index+HALF_SIZE_USER+1,
                                              len(contig_dict["cds_to_keep"]))]
        gc_user_list.extend(gc_chunk)
        target_cds[contig_dict["cds_to_keep"][index]] = gc_chunk
    return gc_user_list, target_cds

def remove_cds_from_dict(list_to_delete, contig_dict):
    for acds in list_to_delete:
        del contig_dict[acds]
    return contig_dict
    
def select_gc_user(target_list, contig_dict):
    contig_infos = ["organism", "strain", "taxon_id", "size"]
    cds_user_list, target_cds = select_user_list(target_list, contig_dict)
    list_to_delete = list(set(contig_dict.keys())-(set(contig_infos)|set(cds_user_list)))
    cds_user_dict = remove_cds_from_dict(list_to_delete, contig_dict)
    cds_user_dict["cds_to_keep"] = cds_user_list
    return cds_user_dict, target_cds

def get_desired_gc(target_list, gc_dict):
    desired_gc = {}
    target_cds = {}
    for acontig in gc_dict:
        desired_gc[acontig], target_cds_contig = select_gc_user(target_list, gc_dict[acontig])
        target_cds.update(target_cds_contig) # rendu possible car le fichier
        # d'input de la boite 2 a été controlé;
        # une target ne peut pas avoir d'informations contradictoires
    return desired_gc, target_cds

def write_json(aFileParsed, filename):
    ''' writes all dictionaries of INSDC file parsed in a json file format
    '''
    with open(filename, "w") as jsonFile:
        json.dump(aFileParsed, jsonFile, indent=4)
    return 0

def main(gcFile, families, targets):
    with open(targets, "r") as target_data:
        targets = target_data.readlines()
    target_list = get_target_list(targets)

    with open(gcFile, "r") as gc_data:
        gc_dict = json.load(gc_data)
    gc_user_dict, target_cds = get_desired_gc(target_list, gc_dict)
    print(target_cds)
    #for cont in gc_user_dict:
    #    print(gc_user_dict[cont].keys())

    #write_json(gc_user_dict, "GC_user.json")

    with open(families, "r") as fam_data:
        fam_dict = json.load(fam_data)
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
