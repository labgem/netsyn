import os
import sys
import json
import subprocess

MICSYNTENY_CP = os.environ.get('MICSYNTENY_CP')
GC_SIZE_USER = 3
HALF_SIZE_USER = int(GC_SIZE_USER/2)
GAP_USER = GC_SIZE_USER
MIN_SYNT_SIZE = 1

def get_target_list(lines):
    target_list = []
    for aline in lines[1:]:
        aline = aline.strip().split("\t")
        target_list.append(aline[0])
    return list(set(target_list))

def find_indices(target_in_contig, cds_list):
    indices = []
    for atarget in target_in_contig:
        indices.append(cds_list.index(atarget))
    return indices

def select_user_list(target_list, contig_dict):
    gc_user_list = []
    cds_target = {}
    target_in_contig = list(set(target_list) & set(contig_dict["cds_to_keep"]))
    indices = find_indices(target_in_contig, contig_dict["cds_to_keep"])
    for index in indices:
        gc_chunk = contig_dict["cds_to_keep"][
            max(index-HALF_SIZE_USER, 0): min(index+HALF_SIZE_USER+1,
                                              len(contig_dict["cds_to_keep"]))]
        gc_user_list.extend(gc_chunk)
        for acds in gc_chunk:
            cds_target[acds] = contig_dict["cds_to_keep"][index]
    return gc_user_list, cds_target

def remove_cds_from_dict(list_to_delete, contig_dict):
    for acds in list_to_delete:
        del contig_dict[acds]
    return contig_dict

def select_gc_user(target_list, contig_dict):
    contig_infos = ["organism", "strain", "taxon_id", "size"]
    cds_user_list, cds_target = select_user_list(target_list, contig_dict)
    list_to_delete = list(set(contig_dict.keys())-(set(contig_infos)|set(cds_user_list)))
    cds_user_dict = remove_cds_from_dict(list_to_delete, contig_dict)
    cds_user_dict["cds_to_keep"] = cds_user_list
    return cds_user_dict, cds_target

def get_desired_gc(target_list, gc_dict):
    desired_gc = {}
    cds_target = {}
    for acontig in gc_dict:
        desired_gc[acontig], cds_target_contig = select_gc_user(target_list, gc_dict[acontig])
        cds_target.update(cds_target_contig) # rendu possible car le fichier
        # d'input de la boite 2 a été controlé;
        # une target ne peut pas avoir d'informations contradictoires
    return desired_gc, cds_target

def get_desired_fam(gc_user_dict, real_fam_dict):
    desired_families = {}
    for afam in real_fam_dict:
        for acont in gc_user_dict:
            desired_families.setdefault(
                afam,
                list(set(real_fam_dict[afam]) & 
                     set(gc_user_dict[acont]["cds_to_keep"])))
            if desired_families[afam] == []:
                del desired_families[afam]
    return desired_families

def write_json(aFileParsed, filename):
    ''' writes all dictionaries of INSDC file parsed in a json file format
    '''
    with open(filename, "w") as jsonFile:
        json.dump(aFileParsed, jsonFile, indent=4)
    return 0

def write_ALL_nodes(cds_target_dict, gc_user_dict):
    INC_pos = 1
    with open("ALL.nodes", "w") as file:
        for acontig in gc_user_dict:
            previous_target = None
            for acds in gc_user_dict[acontig]["cds_to_keep"]:
                target_new = cds_target_dict[acds]
                if not previous_target:
                    previous_target = target_new
                elif target_new != previous_target:
                    INC_pos += GC_SIZE_USER
                    previous_target = target_new
                start_stop_pos = gc_user_dict[acontig][acds]["position"]
                frame = gc_user_dict[acontig][acds]["frame"]
                file.write("\t".join([
                           "@".join([acds, target_new]),
                           str(INC_pos),
                           str(start_stop_pos[0]),
                           str(start_stop_pos[1]),
                           frame,
                           "414"]))
                file.write("\n")
                INC_pos += 1
    return 0

def write_ALL_roles(cds_target_dict, real_fam_dict):
    with open("ALL.roles", "w") as file:
        INC_roles = 1
        for afam in real_fam_dict:
            for member1 in real_fam_dict[afam][:-1]:
                idx = real_fam_dict[afam].index(member1)
                for member2 in real_fam_dict[afam][idx+1:]:
                    target1 = cds_target_dict[member1]
                    target2 = cds_target_dict[member2]
                    file.write(" ".join([
                                str(INC_roles),
                                "@".join([member1, target1]),
                                "@".join([member2, target2])
                                ]))
                    file.write("\n")
                    INC_roles += 1

def run_synteny():
    with open("synteny.log", "w") as file:
        subprocess.run(["java", "-Xmx16G", "-Xms1G", "-classpath", MICSYNTENY_CP, "synteny.synteny", "ALL.nodes", "ALL.nodes", "ALL.roles", "ALL", "-gap", str(GAP_USER), "-size", str(MIN_SYNT_SIZE), ">/dev/null", "||", "exit $?"], stdout=file, stderr=file)
    return 0

def main(gcFile, families, targets):
    with open(targets, "r") as target_data:
        targets = target_data.readlines()
    target_list = get_target_list(targets)
    with open(gcFile, "r") as gc_data:
        gc_dict = json.load(gc_data)
    with open(families, "r") as fam_data:
        fam_dict = json.load(fam_data)

    gc_user_dict, cds_target_dict = get_desired_gc(target_list, gc_dict)
    fam_user_dict = get_desired_fam(gc_user_dict, fam_dict)
    real_fam_dict = {key:value for (key,value) in fam_user_dict.items() if len(value) > 1}

    write_ALL_nodes(cds_target_dict, gc_user_dict)
    write_ALL_roles(cds_target_dict, real_fam_dict)
    
    run_synteny()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
