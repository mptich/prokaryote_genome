# Place a genome (identified by its directory name) into a Taxonomy
# Classification tree, based on the best Common COG correlation
# Directory name is passed in as a parameter

import sys
from taxonomy import *
from collections import defaultdict as DefDict
from common_cogs_method import commonCogsStatDist, buildCogTaxaDict
from shared.algorithms.kendall import calculateWeightedKendall
from shared.pyutils.utils import *
import operator


if len(sys.argv) < 2:
    print("Missing genome directory names")
    sys.exit(-1)

genDirList = sys.argv[1:]
print ("Trying to reclassify taxonomy for %s" % str(genDirList))

cogDict, taxaDict, genCogDict, genTaxaDict = buildCogTaxaDict(3.0, genDirList)

dirCorrList = UtilLoad(GENOME_CORR_LIST())
for t in dirCorrList:
    # Remove all igenomes with correlations < 0.7
    if t[1] < 0.7:
        del cogDict[t[0]]
        del taxaDict[t[0]]

for dir in genCogDict:
    if dir in cogDict:
        del cogDict[dir]
        del taxaDict[dir]

print ("Organisms left in cogDict %d taxaDict %d" %
       (len(cogDict), len(taxaDict)))

taxaTypeSet = set([x.type for x in taxaDict.values()])
print ("taxaTypeSet size %d" % len(taxaTypeSet))

for dir, genCogSet in genCogDict.items():
    print("\nProcessing %s\nCurrent Taxonomy: %s\n" %
          (dir, repr(genTaxaDict[dir])))

    cogDistDict = {}
    for d, cogSet in cogDict.items():
        cogDistDict[d] = commonCogsStatDist(cogSet, genCogSet)
    cogDistList = [cogDistDict[x] for x in sorted(cogDistDict.keys())]

    # Calculate quality of various classifications based on the average
    # COG distance to them
    taxaCogDistDict = DefDict(list)
    for d, taxa in taxaDict.items():
        taxaCogDistDict[taxa.type].append(cogDistDict[d])
    for taxaType, corrList in taxaCogDistDict.items():
        taxaCogDistDict[taxaType] = np.mean(corrList)

    # Now Calculate quality of various classifications based on the
    # correlation between COG distances and taxonomy distances
    taxaCogCorrDict = {}
    for taxaType in taxaTypeSet:
        taxa = Taxa(_name = "Dummy", _type = taxaType)
        taxaDistList = [taxa.distance(taxaDict[x]) for x in
            sorted(taxaDict.keys())]
        corr = calculateWeightedKendall(taxaDistList, cogDistList)
        taxaCogCorrDict[taxaType] = corr

    taxaCogDistList = sorted(taxaCogDistDict.items(),
        key=operator.itemgetter(1))
    taxaCogCorrList = sorted(taxaCogCorrDict.items(),
        key=operator.itemgetter(1), reverse = True)

    print("\nCorrelation List")
    corrLimit = 0.7 * taxaCogCorrList[0][1]
    for taxaType, corr in taxaCogCorrList:
        if corr > corrLimit:
            print("%s: %f" % (repr(taxaType), corr))

    print("\nDistance List")
    distLimit = 1.5 * taxaCogDistList[0][1]
    for taxaType, dist in taxaCogDistList:
        if dist < distLimit:
            print("%s: %f" % (repr(taxaType), dist))





