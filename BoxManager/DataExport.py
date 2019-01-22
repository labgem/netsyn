#!/usr/bin/env python3

##########
# Import #
##########
import common
import logging
#############
# Functions #
#############
def run(logger):
    logger = logging.getLogger('{}.{}'.format(run.__module__, run.__name__))
    logger.info('DataExport running...') ####
    logger.info(common.global_dict) ####