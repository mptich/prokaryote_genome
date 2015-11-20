# This module adds a taxonomy info as taxonomy.json file into all
# genome directories for which we have this info.

import os
import glob
import utils
import csv
import json

# Dictionary that matches genome directory ->
# (taxonomy organism name, name word match count)
writtenTaxonomies = {}

# Taxonomies that could not be matched to any genome
orphanTaxonomies = set()

# Something went wrong
errorCount = 0



dirList = glob.glob(utils.ORGANISMS() + "/*")
dirDict = {}
termDict = {}
for d in dirList:
	_, dir = os.path.split(d)
	dirDict[dir] = set([x.lower() for x in dir.split('_')])
	for t in dirDict[dir]:
		if t == "":
			continue
		dirValList = termDict.get(t, [])
		dirValList.append(dir)
		termDict[t] = dirValList


with open(utils.TAXONOMY(), 'r') as f:
	csvr = csv.reader(f, delimiter = ',')
	for l in csvr:
		if len(l) != 7:
			continue
		name = l[2].lower()
		sname = set(name.split(' '))
		if len(sname) == 0:
			continue
		matched = False
		# Take just one (any) element from this set
		for dir in termDict.get(next(iter(sname)), []):
			if sname.issubset(dirDict[dir]):
				# We got the matching entry in the taxonomy
				# Check if it has been written already
				if dir in writtenTaxonomies:
					# See if it is a better match, if yes - overwrite
					oldName, oldMatchCount = writtenTaxonomies[dir]
					if len(sname) == oldMatchCount:
						# Error - match with the same count
						print(("ERROR: %s is matched again by %s, " +
							"was matched by %s") %
						  	(dir, name, oldName))
						errorCount += 1
						continue

					if (len(sname) < oldMatchCount):
						# This match is worse, just ignore it
						continue

					# We found a better match

				taxonomy = {
					"name": name,
					"id": int(l[0]),
					"domain": l[1],
					"phylum": l[3],
					"class": l[4],
					"order": l[5],
					"family": l[6]
				}
				with open(utils.ORGANISMS() + "/" + dir + "/" +
						utils.TAXONOMY_JSON_FILE(), 'w') as fjson:
					json.dump(taxonomy, fjson)
				writtenTaxonomies[dir] = (name, len(sname))
				matched = True
		if not matched:
			orphanTaxonomies.add(name)

print orphanTaxonomies

print("written %d, orphan %d, errors %d" %
	   (len(writtenTaxonomies), len(orphanTaxonomies), errorCount))






