#!/usr/bin/env python3

##########
# Import #
##########
from netsyn import common

import shutil
import logging
import os
#############
# Functions #
#############


def run(resultsDirectory, analysisNumber, boxes):
    '''
    Updates the completed analysis last number.
    Completes the currently analysis.
    '''
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))

    if not os.path.isdir(resultsDirectory):
        os.mkdir(resultsDirectory)

    allContent = []
    for boxName in boxes[:-1]:
        reportName = common.global_dict['files'][boxName]['report']
        with open(reportName, 'r') as file:
            allContent.append(file.readlines())
    reportName = common.global_dict['files']['EndNetSynAnalysis']['report']
    with open(reportName, 'w') as file:
        for lines in allContent:
            file.write(''.join(lines))

    shutil.copyfile(common.global_dict['files']['DataExport']['graphML'],
                    common.global_dict['files']['EndNetSynAnalysis']['graphML'])
    shutil.copyfile(common.global_dict['files']['DataExport']['html'],
                    common.global_dict['files']['EndNetSynAnalysis']['html'])
    shutil.copytree(common.global_dict['synthesisDataExport'],
                    common.global_dict['synthesisEndNetSynAnalysis'])
    shutil.copyfile(common.global_dict['settingsFileName'],
                    common.global_dict['files']['EndNetSynAnalysis']['settings'])

    with open(common.global_dict['lastAnalysisNumber'], 'w') as file:
        file.write(str(analysisNumber))
    logger.info('Analysis completed, results available in the {} directory'.format(
        resultsDirectory))
