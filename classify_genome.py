#UNFINISHED
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

cogDict, taxaDict, _, _ = buildCogTaxaDict(3.0)
print ("cogDict len %d, taxaDict len %d" % (len(cogDict), len(taxaDict)))

dirCorrDict = UtilLoad(GENOME_CORR_DICT())
print ("dirCorrDict len %d" % len(dirCorrDict))

taxaTypeSet = set([x.type for x in taxaDict.values()])
print ("taxaTypeSet size %d" % len(taxaTypeSet))

print("Building COG distances...")
cogDist = DefDict(dict)
for ordinal, (dir1, cs1) in enumerate(cogDict.iteritems(), start = 1):
    print("\r%d. %s" % (ordinal, dir1)),
    for dir2, cs2 in cogDict.iteritems():
        cogDist[dir1][dir2] = commonCogsStatDist(cs1, cs2)

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





