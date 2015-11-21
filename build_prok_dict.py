# This module builds PROK_GENOME_DICT() file (see utils.py)

import glob
import config
import json
from genome_cls import GenomeError, GenomeJSONEncoder, ProkDna, ProkGenome

# Dictionary that will be saved into the file
prokDict = {}

def newProkGenome(dir):
	prokGenome = ProkGenome(dir)
	fullDir = config.PROKARYOTS_DIR() + dir + '/'
	pttFiles = glob.glob(fullDir + "*.ptt")
	for pttFileName in pttFiles:
		prokDna = ProkDna(pttFileName)
		prokGenome.add(prokDna)
	prokGenome.verify()
	return prokGenome

with open(config.PROKARYOT_DIRS_FILE(), 'r') as fdirs:
	for dir in fdirs:
		dir = dir.strip()
		try:
			prokGenome = newProkGenome(dir)
		except GenomeError as e:
			print("GenomeError: %s" % e)
			continue
		prokDict[dir] = prokGenome

with open(config.PROK_GENOME_DICT(), 'w') as fdict:
	json.dump(prokDict, fdict, cls = GenomeJSONEncoder, sort_keys = True,
			  indent = 4)

