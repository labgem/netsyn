#!/usr/bin/env python3

##########
# Import #
##########
import argparse
import os
import logging
import json
import sys
import math

#############
# Functions #
#############

def find_pos_in_gc(context, prot_ref):
    pos = None
    for idx, gc_member in enumerate(context):
        if gc_member[0] == prot_ref:
            pos = idx+1
            break
    return pos

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

def get_families_per_target(families, cds_info):
    target_famies = {}
    for afam in families:
        for member in families[afam]:
            target = cds_info[str(member[1])]["target"]
            target_famies.setdefault(target, []).append(afam)
    for atarget in target_famies:
        target_famies[atarget] = list(skip_duplicates(target_famies[atarget]))
    return target_famies

def get_syntenies2(families_list, cds_info, families, syntenies):
    for afam in families_list:
        for idx, cdsA in enumerate(families[afam][:-1]):
            for cdsB in families[afam][idx+1:]:
                if not cdsA[0] == cdsB[0]:
                    targA = str(cds_info[str(cdsA[1])]["target"])
                    targB = str(cds_info[str(cdsB[1])]["target"])
                    cdsA_pos = find_pos_in_gc(cds_info[targA]["context"], cdsA[1])
                    cdsB_pos = find_pos_in_gc(cds_info[targB]["context"], cdsB[1])
                    if not (cdsA_pos or cdsB_pos):
                        exit(1)
                    syntenies.setdefault((targA, targB), {}) \
                        .setdefault("syntons", []) \
                        .append((cdsA_pos, cdsB_pos))
                    syntenies[(targA, targB)]["syntons"] = list(skip_duplicates(syntenies[(targA, targB)]["syntons"]))
    return syntenies

def intersect_target_famies(target_famies, cds_info, families):
    syntenies = {}
    targets_list = [target for target in target_famies]
    for idx, targA in enumerate(targets_list[:-1]):
        for targB in targets_list[idx+1:]:
            intersect = list(set(target_famies[targA]).intersection(target_famies[targB]))
            if intersect:
                syntenies = get_syntenies2(intersect, cds_info, families, syntenies)
    return syntenies

# def gc_filter(syntenies, gcUser):
#     for synton in syntenies:
#         for (pos_a, pos_b) in syntenies[synton]:
            

# def synteny_filter(syntenies, gcUser, gap):
#     gc_filter(syntenies, gcUser)
#     #gap_filter()

def run(gcUser, gap, familiesFile, gcFile, contigsFile, TMPDIRECTORY):
    '''
    '''
    BOXNAME = 'SyntenyFinder'
    TMPDIRECTORYPROCESS = '{}/{}'.format(TMPDIRECTORY, BOXNAME)
    if not os.path.isdir(TMPDIRECTORYPROCESS):
        os.mkdir(TMPDIRECTORYPROCESS)
    
    params = {
        "MAX_GC": 5,
        "USER_GC": gcUser
        }
    
    with open(familiesFile, "r") as file:
        families = json.load(file)
    with open(gcFile, "r") as file:
        cds_info = json.load(file)
    with open(contigsFile, "r") as file:
        contigs_info = json.load(file)

    target_famies = get_families_per_target(families, cds_info)
    syntenies = intersect_target_famies(target_famies, cds_info, families)
    print('syntenies:\n{}'.format(syntenies))
#    synteny_filter(syntenies, gcUser, gap)

if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], os.path.abspath("."))


