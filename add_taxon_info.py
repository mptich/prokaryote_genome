# This module adds information about taxonomy for all organisms it can find
# it for

import glob
import csv
from import_proxy import *

# Set of the names of matched taxonomies
matchedTaxonomies = set()

# Taxonomies that could not be matched to any genome
orphanTaxonomies = set()

# Something went wrong
errorCount = 0

# Mapping of name terms -> names, and corresponding ProkDna objects
termDict = {}

with open(config.PROK_CLEAN_GENOME_DICT(), 'r') as fdict:
    masterDict = json.load(fdict, object_hook = UtilJSONDecoderDictToObj)

for d, pds in masterDict.iteritems():
    # Take name from any chromosome; consider what comes before the comma
    pd = pds.getChrom(pds.getChromList()[0])
    name = pd.getName().split(',')[0]
    nameSet = set([x for name.split(' ') if x != ""])
    for term in nameSet:
        o = UtilObject(nameSet = nameSet, name = name, pd = pd)
        l = termDict.get(term, [])
        l.append(o)
        termDict[term] = l

with open(config.TAXONOMY(), 'r') as f:
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
        for obj in termDict.get(next(iter(sname)), []):
            if sname.issubset(obj.nameSet):
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
                with open(config.ORGANISMS() + "/" + dir + "/" +
                        config.TAXONOMY_JSON_FILE(), 'w') as fjson:
                    json.dump(taxonomy, fjson)
                writtenTaxonomies[dir] = (name, len(sname))
                matched = True
        if not matched:
            orphanTaxonomies.add(name)

dirList = glob.glob(config.ORGANISMS() + "/*")
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

print orphanTaxonomies

print("written %d, orphan %d, errors %d" %
       (len(writtenTaxonomies), len(orphanTaxonomies), errorCount))






