# Place a genome (identified by its directory name) into a Taxonomy
# Classification tree, based on the hierarchical algorithm: average
# COG distances for domain, phylum, etc.

import sys
from taxonomy import *
from collections import defaultdict as DefDict
import common_cogs_method as ccm
from shared.algorithms.kendall import calculateWeightedKendall
from shared.pyutils.utils import *
import operator
import math

if len(sys.argv) != 2:
    print "Missing COG distance function name"
    sys.exit(-1)

cogDistFunc = getattr(ccm, sys.argv[1])

cogDict, taxaDict, _, _ = ccm.buildCogTaxaDict(3.0)
print ("cogDict len %d, taxaDict len %d" % (len(cogDict), len(taxaDict)))

dirCorrDict = UtilLoad(GENOME_CORR_DICT())
print ("dirCorrDict len %d" % len(dirCorrDict))

print("Building COG distances...")
cogDist = DefDict(dict)
for ordinal, (dir1, cs1) in enumerate(cogDict.iteritems(), start = 1):
    print("\r%d. %s" % (ordinal, dir1)),
    for dir2, cs2 in cogDict.iteritems():
        cogDist[dir1][dir2] = cogDistFunc(cs1, cs2)

print("\nBuilding average distances for TaxaTypes...")
# Genome dir -> dict of {taxaTypes -> avg COG distance to dir}
dirTaxaTypeDictDict = DefDict(lambda: DefDict(list))
for ordinal, dir1 in enumerate(taxaDict.keys(), start = 1):
    print("\r%d. %s" % (ordinal, dir1)),
    for dir2, taxa in taxaDict.iteritems():
        dirTaxaTypeDictDict[dir1][repr(taxa.type)].append(cogDist[dir1][dir2])

print("\nRebuilding dirTaxaTypeDictDict to get UtilNormDistribs...")
for dir, d in dirTaxaTypeDictDict.iteritems():
    # Find global weighted STD
    std = 0.
    totalLen = 0
    for taxaTypeStr, distList in d.iteritems():
        if len(distList) >= 2:
            val = np.std(distList, ddof = 1.)
            std += val * val * len(distList)
            totalLen += len(distList)
    if totalLen == 0:
        raise ValueError("Cannot calcuate global std for %s" % dir)
    std /= totalLen

    for taxaTypeStr, distList in d.iteritems():
        localStd = np.std(distList)
        localStd *= localStd
        localStd += std
        localStd = math.sqrt(localStd / len(distList))
        d[taxaTypeStr] = UtilNormDistrib(mean=np.mean(distList), std=localStd,
            count=len(distList) + 1)

UtilStore(dirTaxaTypeDictDict, TAXA_TYPE_COG_DIST_DICT())

# Now that we have build dictionary of average distances, let's translate
# it into the best reclassified TaxaTypes

# First, build dictionary dir -> TaxaTypeTree
print ("\nBuilding TaxaType trees...")
taxaTypeTree = TaxaTypeTree(set([x.type for x in taxaDict.values()]))
UtilStore(taxaTypeTree, TAXA_TYPE_TREE())

taxaDistCntDict = UtilLoad(GENOME_TAX_DIST_CNT_DICT())

# Now find optimal reclassified TaxaTypes, and dump them into a file
print("Build reclassification...")
reclassObjList = []
dumpDirNodeCostDict = {}
for dir in dirTaxaTypeDictDict.keys():
    nodeCostDict = taxaTypeTree.bldCostDict(dirTaxaTypeDictDict[dir])
    dumpDirNodeCostDict[dir] = taxaTypeTree.utilJsonDump(
        nodeAttribDict = nodeCostDict)
    taxaType, cost = taxaTypeTree.optimal(nodeCostDict)
    dist = taxaDict[dir].type.distance(taxaType)
    reclassObjList.append(UtilObject(dir = dir, cogCorr = dirCorrDict[dir],
        oldClassif = taxaDict[dir].type, newClassif = taxaType,
        taxaDist = dist, cogDist=cost, taxaDistCnts = taxaDistCntDict[dir]))

UtilStore(dumpDirNodeCostDict, DIR_NODE_COST_DICT())

reclassObjList = sorted(reclassObjList, key = lambda x: x.cogCorr)

UtilStore([x for x in reclassObjList if x.taxaDist > 0],
          HIER_RECLASSIFIED_LIST())

distList = [0] * (TaxaType.maxDistance() + 1)
for obj in reclassObjList:
    distList[obj.taxaDist] += 1
print("Out of %d genomes, reclassification dist distribution %s" %
      (len(dirTaxaTypeDictDict), repr(distList)))

# Calculate Kendal correlation between taxaDist and cogCorr
corr = calculateWeightedKendall([x.taxaDist for x in reclassObjList],
    [x.cogCorr for x in reclassObjList])
print "taxaDist / cogCorr correlation", corr