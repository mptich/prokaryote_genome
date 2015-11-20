# This module builds GENOME_DICT_JSON file.
# This file maps directory -> {strain1: <DNA in chromosomes>, strain2: <DNA in chromosomes>, ...,
# plasmids: [list of plasmid DNA objects], phages: [list of phage DNA objects] }
# [list of plasmid genomes]}
# Each genome is described by: {ptt: PTT file name, name: genome's name}

import os
import glob
import utils
import csv
import json
import re


# Dictionary that will be saved into the file
genDict = {}

# Pattern of plasmit string
plasmidPat = re.compile(r'.*plasmid.*')

def processPttFile(dir, pttName):
	global genDict
	d = genDict.get(dir, {"plasmids":[]})

	# Get just the first line from the ptt file
	with open(pttName, 'r') as fptt:
		for l in fptt:
			# Get it up to a comma
			genName = l.split(',')[0]
			break

	genEntry = {"ptt": os.path.split(pttName)[1], "name": genName}
	if plasmidPat.match(genName):
		# This is a plasmid DNA
		d["plasmids"].append(genEntry)
	else:
		# This is the main DNA
		if "genome" in d:
			# Second genome file
			print("ERROR: %s got genomes %s and %s" %
				  (dir, d["genome"]["name"], genName))
		else:
			d["genome"] = genEntry
	genDict[dir] = d


with open(utils.GENOME_DIRS_FILE(), 'r') as fdirs:
	for fullDir in fdirs:
		# Remove new line
		fullDir = fullDir.rstrip()
		pttList = glob.glob(fullDir + "/*.ptt")
		for pttName in pttList:
			processPttFile(os.path.split(fullDir)[1], pttName)

with open(utils.GENOME_DICT_JSON(), 'w') as fjson:
	json.dump(genDict, fjson)
