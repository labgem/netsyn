#!/usr/bin/env python3

##########
# Import #
##########
import shutil
import common
import logging
#############
# Functions #
#############
def run(resultsDirectory, analysisNumber):
    '''
    Updates the completed analysis last number.
    Completes the currently analysis.
    '''
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    shutil.copyfile(common.global_dict['settingsFileName'], common.global_dict['files']['EndNetSynAnalysis']['settings'])
    with open(common.global_dict['lastAnalysisNumber'], 'w') as file:
        file.write(str(analysisNumber))
    logger.info('Analysis completed, results available in the {} directory'.format(resultsDirectory))
