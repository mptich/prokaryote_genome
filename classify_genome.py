# Place a genome (identified by its directory name) into a Taxonomy
# Classification tree, based on the hierarchical algorithm: average
# COG distances for domain, phylum, etc.

import sys
from taxonomy import *
from collections import defaultdict as DefDict
import common_cogs_method as commonCogsMethod
from shared.algorithms.kendall import calculateWeightedKendall
from shared.pyutils.utils import *
#from shared.pyutils.UtilNormDistrib import *
import config
import operator
import math
import itertools

CutOffDiff = 1.0
CutOffBestFit = 3.0


_, _, taxaDict, _ = \
    commonCogsMethod.buildCogTaxaDict(noWeights = True)
print ("taxaDict len %d" % len(taxaDict))

print("Reading COG distances...")
cogDist = UtilLoad(COG_DIST_DICT())

# Build a tree of TaxaTypes
taxaTypeTree = TaxaTypeTree(taxaDict)

# Set of all Taxa types on all levels
allTaxaTypes = taxaTypeTree.getAllTypesSet()
print("Length of allTaxaTypes %d" % len(allTaxaTypes))

# Build a dictionary: [dir][taxaType] -> UtilObject(mean, std,
# isAncest, distList), where
# mean - mean distance between this dir and all [other] dirs in this taxaType
# std - standard deviation, if applicable
# distList - list of all the distances
print("Building taxaTypeDistDict...")
taxaTypeDistDict = DefDict(dict)
taxaTypeAncDistDict = DefDict(lambda: [None] * (TaxaType.hierarchySize()+1))
for ind, dir in enumerate(taxaDict, start = 1):
    print("\r%u. %s" % (ind, dir)),
    d = taxaTypeDistDict[dir]
    listAnc = taxaTypeAncDistDict[dir]
    cogDistForDir = cogDist[dir]
    for currType in allTaxaTypes:
        dirs = taxaTypeTree.getDirSet(currType)
        if dir in dirs:
            ancestorList = True
            dirs.remove(dir)
        else:
            ancestorList = False
        distList = [cogDistForDir[x] for x in dirs]
        meanVal = 0.0
        stdVal = 0.0
        if len(distList) > 0:
            meanVal = np.mean(distList)
        if len(distList) > 1:
            stdVal = np.std(distList, ddof = 1.0)
        if ancestorList:
            listAnc[currType.depth()] = UtilObject(mean=meanVal, std=stdVal,
                distList=distList)
        else:
            d[currType] = UtilObject(mean=meanVal, std=stdVal,
                distList=distList)
print

# Build lists of distances for each level
print("Build globDistList...")
globDistList = []
for i in range(TaxaType.hierarchySize() + 1):
    globDistList.append([])
for listAnc in taxaTypeAncDistDict.values():
    for ind, obj in enumerate(listAnc):
        if obj is None:
            continue
        globDistList[ind].extend(obj.distList)

# Build list of UtilObject(mean, std, count)
globStdList = []
for l in globDistList:
    UtilDrawHistogram(l, show = False)
    if len(l) >= 2:
        std = std=np.std(l, ddof=1.0)
    else:
        std = None
    globStdList.append(std)
UtilDrawHistogram(show = True)
print globStdList

bestFitHistogram = []
reclassList = []
print("RECLASSIFICATIONS...")
for ind, (dir, taxa) in enumerate(taxaDict.iteritems(), start=1):
    typeOrig = taxa.type
    distObjAnc = taxaTypeAncDistDict[dir]
    taxaTypeToObj = taxaTypeDistDict[dir]

    bestFit = 0.0
    bestFitType = None
    bestFitComparedTaxons = None
    for taxaOther in taxaDict.values():
        typeOther = taxaOther.type
        commonAncestor = typeOrig.commonAncestor(typeOther)
        commonDepth = commonAncestor.depth()
        if (commonAncestor == typeOrig) or (commonAncestor == typeOther):
            continue
        distObjOther = [None] * (TaxaType.hierarchySize() + 1)
        currType = typeOther
        while currType.depth() > commonDepth:
            distObjOther[currType.depth()] = taxaTypeToObj[currType]
            currType = currType.parent()
        sum = 0.0
        worstDiff = 1000000. # Very large number
        compCount = 0
        comparedTaxons = []
        prevDistAnc = None
        prevDistOther = None
        for i in range(TaxaType.hierarchySize(), commonDepth, -1):
            objAnc = distObjAnc[i]
            objOther = distObjOther[i]
            calcAnc = bool(objAnc) and bool(objAnc.distList)
            calcOther = bool(objOther) and bool(objOther.distList)
            if calcAnc:
                prevDistAnc = objAnc.mean
            if calcOther:
                prevDistOther = objOther.mean
            if (calcAnc or calcOther) and \
                (prevDistAnc is not None) and (prevDistOther is not None):
                diff = (prevDistAnc - prevDistOther) / globStdList[i]
                if diff < CutOffDiff:
                    sum = 0.0
                    break
                sum += diff
                comparedTaxons.append((TaxaType.hierarchy()[i-1], diff))
        if sum > bestFit:
            bestFit = sum
            bestFitType = typeOther
            bestFitComparedTaxons = comparedTaxons

    print("\r%u. %s bestFit %f" % (ind, dir, bestFit)),
    bestFitHistogram.append(bestFit)
    if bestFit >= CutOffBestFit:
        s = ("\n%s\nOriginal: %s\nReclassified: %s\nTaxonomy distance: %d" +\
            "\nSigmas: %f\nCompared taxons: %s\nSigmas per compare: %f\n")%\
            (dir, repr(typeOrig), repr(bestFitType),
            typeOrig.distance(bestFitType), bestFit,
            ", ".join([":".join((str(y) for y in x)) for x in \
            bestFitComparedTaxons]),
            bestFit/len(bestFitComparedTaxons))
        print s
        reclassList.append((bestFit, s))

UtilDrawHistogram(bestFitHistogram, show=True)

reclassList = sorted(reclassList, reverse = True)

with open(config.WORK_FILES_DIR() + "Reclassify.txt", "w") as f:
    for t in reclassList:
        f.write(t[1])
