# This module builds PROK_GENOME_DICT() file (see utils.py)

import glob
import utils
import json
from genome_cls import GenomeError, ProkDna, ProkGenome

# Dictionary that will be saved into the file
prokDict = {}

def newProkGenome(dir):
	prokGenome = ProkGenome(dir)
	fullDir = utils.PROKARYOTS_DIR() + dir + '/'
	pttFiles = glob.glob(fullDir + "*.ptt")
	for pttFileName in pttFiles:
		prokDna = ProkDna(pttFileName)
		prokGenome.add(prokDna)
	prokGenome.verify()
	return prokGenome

with open(utils.PROKARYOT_DIRS_FILE(), 'r') as fdirs:
	for dir in fdirs:
		dir = dir.strip()
		try:
			prokGenome = newProkGenome(dir)
		except GenomeError as e:
			print("GenomeError: %s" % e)
			continue
		prokDict[dir] = prokGenome

with open(utils.PROK_GENOME_DICT(), 'w') as fdict:
	json.dump(prokDict, fdict, sort_keys = True, indent = 4)

