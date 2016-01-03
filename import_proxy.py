# This module makes sure that all modules imported from dynamic locations
# (like shared) are loaded only once
# Also it abstracts imp module, which has been deprecated in Python 3.4

import config
import os
import imp
import json

def IMPORT_BY_PATH(path):
    modName = os.path.basename(path).split('.')[0]
    return imp.load_source(modName, path)

IMPORT_BY_PATH(config.SHARED_PROG_DIR() + "pyutils/utils.py")
from utils import *

IMPORT_BY_PATH(config.SHARED_PROG_DIR() + "algorithms/kendall.py")
from kendall import calculateWeighetedKendall

from genome_cls import ProkDna, ProkDnaSet, ProkGenome, CogInst, Cog

