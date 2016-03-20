# Place a genome (identified by its directory name) into a Taxonomy
# Classification tree, based on the best Common COG correlation
# Directory name is passed in as a parameter

import sys
from taxonomy import *
from collections import defaultdict as DefDict
from common_cogs_method import commonCogsStatDist, buildCogTaxaDict
from shared.algorithms.kendall import calculateWeightedKendall
from shared.pyutils.utils import *


if len(sys.argv) < 2:
    print("Missing genome directory names")
    sys.exit(-1)

genDirList = sys.argv[1:]
print ("Placing genomes under directories %s" % str(genDirList))

cogDict, taxaDict, genDict = buildCogTaxaDict(3.0, genDirList)

dirCorrList = UtilLoad(GENOME_CORR_LIST())
for t in dirCorrList:
    # Remove all igenomes with correlations < 0.7
    if t[1] < 0.7:
        del cogDict[t[0]]
        del taxaDict[t[0]]

for dir in genDict:
    if dir in cogDict:
        del cogDict[dir]
        del taxaDict[dir]

print ("Organisms left in cogDict %d taxaDict %d" %
       (len(cogDict), len(taxaDict)))

taxaTypeSet = set([x.type for x in taxaDict.values()])
print ("taxaTypeSet size %d" % len(taxaTypeSet))

for dir, genCogSet in genDict.items():
    print("\nProcessing %s" % dir)

    cogDistList = [commonCogsStatDist(cogDict[d], genCogSet) for d in
        sorted(cogDict.keys())]

    bestCorr = -1.1
    for taxaType in taxaTypeSet:
        taxa = Taxa(_name = "Dummy", _type = taxaType)
        taxaDistList = [taxa.distance(taxaDict[x]) for x in
            sorted(taxaDict.keys())]
        corr = calculateWeightedKendall(taxaDistList, cogDistList)
        if corr > bestCorr:
            bestCorr = corr
            bestTaxaType = taxaType
            print("bestCorr %f taxaType %s" % (bestCorr, str(taxaType)))

