##########
# Import #
##########
import __main__ as namespace

#############
# Functions #
#############
def unPetitBonjourPourredonnerLeMoral(msg=None):
    '''
    Fonction de desespoire pour les developpeur en quete de debogage...
    '''
    print("=D I'M HAPPY!!!")
    if msg:
        print(msg)


def constantsInitialiszation(projectName, inputFile,):
    '''
    Initialization of constants.
    Calling from netsyn or BoxManager modules.
    '''
    global_dict['workingDirectory'] = projectName
    global_dict['tmpDirectory'] = '{}/TMP'.format(projectName) #TMPDIRECTORY
    global_dict['settingsFileName'] = '{}/{}'.format(projectName, '.lastSettings.yml') #SETTINGSFILENAME
    global_dict['reportFileName'] = '{}/{}'.format(projectName, '.report.yml') #REPORTFILENAME
    global_dict['inputFileSaved'] = '{}/{}'.format(projectName, inputFile) #INPUTLIST

def filesNameInitialiszation(resultsDirectory, projectName, analysisNumber):
    global_dict['files'] = {
        global_dict['boxName']['GetINSDCFiles'] : {
            'inputClusteringStep' : '{}/{}/inputClusteringIntoFamiliesStep.tsv'.format(global_dict['tmpDirectory'], global_dict['boxName']['GetINSDCFiles'])
        },
        global_dict['boxName']['ClusteringIntoFamilies'] : {
            'faa' : '{}/{}/MMseqs2_run.faa'.format(global_dict['tmpDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'contigs' : '{}/{}/contigs.pickle'.format(global_dict['tmpDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'genomicContexts' : '{}/{}/genomicContexts.pickle'.format(global_dict['tmpDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'lineage' : '{}/{}/taxonomyLineage.pickle'.format(global_dict['tmpDirectory'], global_dict['boxName']['ClusteringIntoFamilies']),
            'targets' : '{}/{}/targets_list'.format(global_dict['tmpDirectory'], global_dict['boxName']['ClusteringIntoFamilies'])
        },
        global_dict['boxName']['SyntenyFinder'] : {
            'clustering' : '{}/{}/Clustering.json'.format(global_dict['tmpDirectory'], global_dict['boxName']['SyntenyFinder'])
        },
        global_dict['boxName']['DataExport'] : {
            'graphML' : '{}/{}_Results_{}.graphML'.format(resultsDirectory, projectName, analysisNumber),
            'html' : '{}/{}_Results_{}.html'.format(resultsDirectory, projectName, analysisNumber),
            'settings' : '{}/{}_Settings_{}.yaml'.format(resultsDirectory, projectName, analysisNumber)
        }
    }


#########################
# Constantes definition #
#########################
global_dict = {
    'defaultValue': 'NA',
    'maxGCSize': 11, #MAXGCSIZE
    'boxName' : {
        'GetINSDCFiles' : 'GetINSDCFiles',
        'ClusteringIntoFamilies' :'ClusteringIntoFamilies',
        'SyntenyFinder' : 'SyntenyFinder',
        'DataExport' : 'DataExport'
    },
    'inputIIheaders' : [
        'UniProt_AC',
        'protein_AC',
        'protein_AC_field',
        'nucleic_AC',
        'nucleic_File_Format',
        'nucleic_File_Name'
    ]
}


##############################
# Add variables to namespace #
##############################
for variable,value in global_dict.items():
    setattr(namespace, variable, value)