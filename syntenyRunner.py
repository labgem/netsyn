import os
import sys
import json
import subprocess
from py4j.java_gateway import JavaGateway

MICSYNTENY_CP = "/".join([sys.path[0], "synteny"])
#MICSYNTENY_CP = "/".join([sys.path[0], "jaSynt"]) # equivalent cmd : os.path.dirname(os.path.abspath(__file__))
#MICSYNTENY_CP = os.getenv('MICSYNTENY_CP')
#print(MICSYNTENY_CP)
GC_SIZE_USER = 5
HALF_SIZE_USER = int(GC_SIZE_USER/2)
GAP_USER = 3
MIN_SYNTENY_SIZE = 1
FILTERING = "ON"
analysis_id = "414"

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
        cds_target[contig_dict["cds_to_keep"][index]] = gc_chunk
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

def set_relative_position(cds_target_dict, gc_user_dict):
    target_cds_posRel = {}
    INC_pos = 1
    all_targets_ordered = sorted(cds_target_dict.keys())
    contigs_size = [len(gc_user_dict[acontig]["cds_to_keep"])
                    for acontig in gc_user_dict]
    all_contig = dict(zip(gc_user_dict.keys(), contigs_size))

    for atarget in all_targets_ordered:
        for acontig in all_contig:
            if atarget in gc_user_dict[acontig]["cds_to_keep"]:
                target_cds_posRel[atarget] = []
                for acds in cds_target_dict[atarget]:
                    target_cds_posRel[atarget].extend(
                        list(zip((acds,), (str(INC_pos),))))
                    INC_pos += 1
            if all_contig[acontig] > 1:
                all_contig[acontig] = all_contig[acontig]-1
            else:
                del all_contig[acontig]
        INC_pos += GC_SIZE_USER
    return target_cds_posRel


def write_ALL_nodes(gc_user_dict, posRel_dict):
    with open("ALL.nodes", "w") as file:
        for atarget in sorted(posRel_dict.keys()):
            for acontig in gc_user_dict:
                if atarget in gc_user_dict[acontig]["cds_to_keep"]:
                    for acds, aposRel in posRel_dict[atarget]:
                        start, stop = gc_user_dict[acontig][acds]["position"]
                        frame = gc_user_dict[acontig][acds]["frame"]
                        file.write("\t".join([
                            "@".join([acds, atarget]),
                            str(aposRel),
                            str(start),
                            str(stop),
                            frame,
                            analysis_id]))
                        file.write("\n")
    return 0

def store_roles(real_fam_dict, cds_target_dict):
    roles_storage = []
    INC_roles = 1
    for afam in real_fam_dict:
        for member1 in real_fam_dict[afam][:-1]:
            idx = real_fam_dict[afam].index(member1)
            for member2 in real_fam_dict[afam][idx+1:]:
                target1 = [atarget for atarget in cds_target_dict if member1
                           in cds_target_dict[atarget]]
                target2 = [atarget for atarget in cds_target_dict if member2
                           in cds_target_dict[atarget]]
                for atarget1 in target1:
                    for atarget2 in target2:
                        roles_storage.append((INC_roles, member1, atarget1,
                                              member2, atarget2))
                        INC_roles += 1
    return roles_storage

def write_ALL_roles(roles_storage):
    with open("ALL.roles", "w") as file:
        for arole in roles_storage:
            file.write(" ".join([
                str(arole[0]),
                "@".join([arole[1], arole[2]]),
                "@".join([arole[3], arole[4]])
                ]))
            file.write("\n")
    return 0

def run_synteny():
    LOCAL_FILE_PATH = os.path.abspath(os.path.dirname(__file__))
    SYNTENYJAR = "/".join([LOCAL_FILE_PATH, "jar/synteny.jar"])
    with open("synteny.log", "w") as file:
        synteny_run = subprocess.run(["java", "-Xmx16G", "-Xms1G", "-classpath", SYNTENYJAR,
                                      "synteny.synteny", "ALL.nodes", "ALL.nodes", "ALL.roles",
                                      "ALL", "-gap", str(GAP_USER), "-size",
                                      str(MIN_SYNTENY_SIZE)],
                                     stdout=file, stderr=file, check=True)
        print('exit {}: {}'.format("synteny_run", synteny_run.returncode))
    return 0

    # with open("synteny.log", "w") as file:
    #     subprocess.run(["java", "-Xmx16G", "-Xms1G", "-classpath", MICSYNTENY_CP,
    #                     "synteny.synteny", "ALL.nodes", "ALL.nodes", "ALL.roles",
    #                     "ALL", "-gap", str(GAP_USER), "-size",
    #                     str(MIN_SYNTENY_SIZE), ">/dev/null", "||", "exit $?"],
    #                    stdout=file, stderr=file)
    # return 0

def get_cpd_role_synton(lines):
    role_synton_cpd = {}
    role_synton_cpd["roles"] = {}
    role_synton_cpd["syntons"] = {}
    for aline in lines:
        aline = aline.strip().split("$")
        role, synton = aline
        role_synton_cpd["roles"][role] = synton
        role_synton_cpd["syntons"].setdefault(synton, []).append(role)
    print("role_synton_cpd :", role_synton_cpd)
    return role_synton_cpd

def synteny_filtering_preparation(posRel_dict, roles_storage, role_synton_cpd):
    pos_in_synteny = []
    for asynteny in role_synton_cpd["syntons"]:
        for aroleID in role_synton_cpd["syntons"][asynteny]:
            role = [arole for arole in roles_storage if str(arole[0]) == aroleID]
            cds1, target1, cds2, target2 = role[0][1:]
            for atarget in posRel_dict:
                if atarget == target1:
                    for cds, position in posRel_dict[atarget]:
                        if cds == cds1:
                            pos_cds1 = position
                        if cds == target1:
                            pos_target1 = position
                if atarget == target2:
                    for cds, position in posRel_dict[atarget]:
                        if cds == cds2:
                            pos_cds2 = position
                        if cds == target2:
                            pos_target2 = position
            pos_in_synteny.append((asynteny, pos_cds1, pos_target1, pos_cds2, pos_target2))
    return pos_in_synteny

def min_max_synteny(pos_in_synteny):
    synteny_description = {}
    for asynteny in pos_in_synteny:
        synt_id = asynteny[0]
        pos_gene_1 = asynteny[1]
        pos_target_gene_1 = asynteny[2]
        pos_gene_2 = asynteny[3]
        pos_target_gene_2 = asynteny[4]

        if not synt_id in synteny_description:
            synteny_description[synt_id] = {'minG1' : pos_gene_1,
                                            'maxG1' : pos_gene_1,
                                            'targetG1' : pos_target_gene_1,
                                            'minG2' : pos_gene_2,
                                            'maxG2' : pos_gene_2,
                                            'targetG2' : pos_target_gene_2}
        else:
            # Genome 1 in synteny
            if pos_gene_1 < synteny_description[synt_id]['minG1']:
                synteny_description[synt_id]['minG1'] = pos_gene_1
            elif pos_gene_1 > synteny_description[synt_id]['maxG1']:
                synteny_description[synt_id]['maxG1'] = pos_gene_1
            # Genome 2 in synteny
            if pos_gene_2 < synteny_description[synt_id]['minG2']:
                synteny_description[synt_id]['minG2'] = pos_gene_2
            elif pos_gene_2 > synteny_description[synt_id]['maxG2']:
                synteny_description[synt_id]['maxG2'] = pos_gene_2
    return synteny_description

def genome_conserved(start, end, target, gap):
    if target >= start and target <= end:
        conserved = True                        # .....start....target....end.....
    elif target < start:
        # pas besoin du if-else:
        # conserved = (start-target-1) <= gap
        if (start-target-1) <= gap:
            conserved = True                    # .....target.start....end.....
        else:
            conserved = False                   # .....target...........start....end.....
    elif target > end:
        # idem
        # conserved = (target-end-1) <= gap
        if (target-end-1) <= gap:
            conserved = True                    # .....start.....end.target.....
        else:
            conserved = False                   # .....start.....end...........target.....
    return conserved

def check_target_ourside_synteny(mini, maxi, target):
    if target < mini:
        mini = target
    elif target > maxi:
        maxi = target
    return mini, maxi

def synteny_purification(synteny_description, gap, filtering):
    synteny_conserved = []
    synteny_not_conserved = []
    for synt_id in synteny_description:
        minG1 = synteny_description[synt_id]['minG1']
        maxG1 = synteny_description[synt_id]['maxG1']
        targetG1 = synteny_description[synt_id]['targetG1']
        minG2 = synteny_description[synt_id]['minG2']
        maxG2 = synteny_description[synt_id]['maxG2']
        targetG2 = synteny_description[synt_id]['targetG2']

        if filtering == 'on':
            g1 = genome_conserved(minG1, maxG1, targetG1, gap)
            g2 = genome_conserved(minG2, maxG2, targetG2, gap)
        else:
            g1 = True
            g2 = True

        minG1, maxG1 = check_target_ourside_synteny(minG1, maxG1, targetG1)
        minG2, maxG2 = check_target_ourside_synteny(minG2, maxG2, targetG2)

        if g1 and g2:
            synteny_conserved.append([synt_id, minG1, maxG1, minG2, maxG2])
        else:
            synteny_not_conserved.append([synt_id, minG1, maxG1, minG2, maxG2])

    return synteny_conserved, synteny_not_conserved

def list_writing(list, output_name, separator):
    with open(output_name, 'w') as file:
        for line in list:
            line = separator.join(str(e) for e in line)
            file.write(line+'\n')

def outputs_writing(synteny_conserved, synteny_not_conserved, separator):
    list_writing(synteny_conserved, 'synteny_conserved.filtered', separator)
    list_writing(synteny_not_conserved, 'synteny_not_conserved.filtered', separator)

def filter_synteny(positions_in_synteny):
    # Get min and max relative positions for Genome 1 and 2 in synteny
    synteny_description = min_max_synteny(positions_in_synteny)
    # Check filter
    synteny_conserved, synteny_not_conserved = synteny_purification(synteny_description, GAP_USER, FILTERING)
    # Outputs writing
    separator = "$"
    outputs_writing(synteny_conserved, synteny_not_conserved, separator)
    return synteny_conserved

def get_syntons_info(synton_db):
    synton_desc = {}
    for aline in synton_db:
        aline = aline.strip().split("$")
        synt_id = aline[0]
        analysis_id = aline[1]
        startG1 = aline[3]
        stopG1 = aline[4]
        startG2 = aline[5]
        stopG2 = aline[6]
        nbgeneG1 = aline[7]
        nbgeneG2 = aline[8]
        synton_desc[synt_id] = [
            synt_id,
            analysis_id,
            startG1,
            stopG1,
            startG2,
            stopG2,
            nbgeneG1,
            nbgeneG2
            ]
    return synton_desc

def score_calculation(syntenies_list, synton_desc):
    syntenies_scores = {}
    for asynt_id in synton_desc:
        prot_in_synt = int(synton_desc[asynt_id][6]) + int(synton_desc[asynt_id][7])
        startG1, stopG1, startG2, stopG2 = [pos for pos in syntenies_list if pos[0] == asynt_id][0][1:]
        synteny_rangeG1 = (int(stopG1)-int(startG1))+1
        synteny_rangeG2 = (int(stopG2)-int(startG2))+1
        mean_score = prot_in_synt/2
        weighted_score = prot_in_synt / (synteny_rangeG1 + synteny_rangeG2)
        syntenies_scores.setdefault(asynt_id, {}).setdefault("mean", mean_score)
        syntenies_scores[asynt_id]["weighted"] = weighted_score
    return syntenies_scores

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
    real_fam_dict = {key:value for (key, value) in fam_user_dict.items() if len(value) > 1}
    posRel_dict = set_relative_position(cds_target_dict, gc_user_dict)
    roles_storage = store_roles(real_fam_dict, cds_target_dict)
    
    write_ALL_nodes(gc_user_dict, posRel_dict)
    write_ALL_roles(roles_storage)
    run_synteny()

    with open("ALL.CPD.DB", "r") as cpd_id_roles_synton:
        id_roles_synton = cpd_id_roles_synton.readlines()
    role_synton_cpd = get_cpd_role_synton(id_roles_synton)

    positions_in_synteny = synteny_filtering_preparation(posRel_dict, roles_storage, role_synton_cpd)
    conserved_synteny = filter_synteny(positions_in_synteny)

    with open("ALL.SYNTON.DB", "r") as synton_db:
        short_synton_db = synton_db.readlines()
    synton_desc = get_syntons_info(short_synton_db)

    syntenies_scores = score_calculation(conserved_synteny, synton_desc)
    #write_scores(mean_score, weighted_score,) # pas utile d'éditer le fichier
                                               # ALL.SYNTON.DB car pas de
                                               # nouvelle entrée dans la base
                                               # de données
    print("\n", synton_desc)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
    #print("\t".join([var for var in globals() if var[0] != "_"]))

