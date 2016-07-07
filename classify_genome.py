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


if len(sys.argv) != 2:
    print "Missing COG distance function name"
    sys.exit(-1)

cogDistFuncName = sys.argv[1]

_, _, taxaDict =\
    commonCogsMethod.buildCogTaxaDict()
print ("taxaDict len %d" % len(taxaDict))

print("Reading COG distances...")
cogDist = UtilLoad(COG_DIST_DICT(cogDistFuncName))

# Build dictionaries similar to taxaDict, but for all variants of
# depths of the taxonomy
print("Building taxaTypeDictList...")
taxaTypeDictList = [DefDict(lambda: None) for i in range(
    TaxaType.maxDistance()+1)]
for dir, taxa in taxaDict.iteritems():
    taxaType = taxa.type
    depth = taxaType.depth()
    while True:
        taxaTypeDictList[depth][dir] = taxaType
        if depth == 0:
            break
        taxaType = taxaType.parent()
        depth = taxaType.depth()

# Now build dictionaries of taxa type -> set of dirs, for each depth level
print("Building taxaTypeRevDictList...")
taxaTypeRevDictList = [DefDict(set) for i in range(TaxaType.maxDistance()+1)]
for ind, ttd in enumerate(taxaTypeDictList):
    d = taxaTypeRevDictList[ind]
    for dir, taxaType in ttd.iteritems():
        d[taxaType].add(dir)

# Calculate global distribution of distances, on each level
print("Distribution of distances...")
distDistrDirList = [DefDict(lambda: DefDict(list)) for i in \
    range(TaxaType.maxDistance()+1)]
for ind, (dir1, d) in enumerate(cogDist.iteritems(), start=1):
    print("\r%u. %s" % (ind, dir1)),
    for dir2, dist in d.iteritems():
        if dir1 == dir2:
            continue
        # Find common parent
        for depth in reversed(range(TaxaType.maxDistance()+1)):
            taxaType1 = taxaTypeDictList[depth][dir1]
            taxaType2 = taxaTypeDictList[depth][dir2]
            if taxaType1 and taxaType2 and (taxaType1 == taxaType2):
                distDistrDirList[depth][taxaType1][dir1].append(dist)
                break
print

print("Building distDistrList...")
distDistrList = [DefDict(list) for i in range(TaxaType.maxDistance()+1)]
for ind, dd in enumerate(distDistrDirList):
    print("Size of dict at depth %d in distDistrDirList is %d" %
          (ind, len(dd)))
    for taxaType, dird in dd.iteritems():
        for dir, distList in dird.iteritems():
            distDistrList[ind][taxaType].append(np.mean(distList))

print("Histogram of distances...")
for depth in range(TaxaType.maxDistance()+1):
    inputList = list(itertools.chain.from_iterable(
        distDistrList[depth].values()))
    if depth != 0:
        UtilDrawHistogram(inputList = inputList, show = False)
UtilDrawHistogram(show=True)

# Calculate global STD. STD is calculated for each level of the Taxonomy
# tree separately
print("Calculating globalStd...")
globStd = []
for depth in range(TaxaType.maxDistance()+1):
    sumInst = 0
    sumStd = 0.0
    for taxaType, distList in distDistrList[depth].iteritems():
        if len(distList) < 3:
            continue
        sumInst += len(distList)
        std = np.std(distList, ddof = 1)
        sumStd += len(distList) * std * std
    globStd.append(math.sqrt(sumStd / sumInst))
print("globStd: %s" % globStd)

# Build genome reclassifications
tempSameCount = 0
tempGoodSigmaCount = 0
tempList = []
tempObjList = []
print("Build reclassification...")
for ind, (dir, taxa) in enumerate(taxaDict.iteritems(), start=1):
    print("\r%u. %s" % (ind, dir)),
    oldType = taxa.type
    origType = oldType
    depth = oldType.depth()
    while len(taxaTypeRevDictList[depth][oldType]) < 2:
        oldType = oldType.parent()
        depth = oldType.depth()
    # Find old distance
    distList = [cogDist[dir][dir1] for dir1 in \
        taxaTypeRevDictList[depth][oldType] if dir != dir1]
    oldDist = np.mean(distList)
    # Find better type
    betterDist = oldDist
    betterType = None
    for taxaType, dirSet in taxaTypeRevDictList[depth].iteritems():
        if taxaType == oldType:
            continue
        assert(dir not in dirSet)
        distList = [cogDist[dir][dir1] for dir1 in dirSet]
        newDist = np.mean(distList)
        if newDist < betterDist:
            betterDist = newDist
            betterType = taxaType

    if not betterType:
        tempSameCount += 1
    else:
        # Reclassified !
        assert(oldType.depth() == betterType.depth())
        depth = betterType.depth()
        sigmas = (oldDist - betterDist) / globStd[depth]
        if sigmas >= 3.0:
            tempGoodSigmaCount += 1
        desc = "%s\n" \
            "Taxonomy distance between old and new positions: %d\n" \
            "Sigmas: %f\n" \
            "Old classification: %s\n" \
            "New classification: %s\n" \
            "Old COG distance: %f\n" \
            "New COG distance: %f\n" \
            "STD: %f\n\n" % \
            (dir, origType.distance(betterType), sigmas,
            repr(origType), repr(betterType), oldDist,
            betterDist, globStd[depth])
        tempList.append((sigmas, oldDist - betterDist, desc))
        obj = UtilObject(dir = dir,
            taxaDist = origType.distance(betterType), sigmas=sigmas, \
            oldClass = origType, newClass=betterType, oldCogDist = oldDist, \
            newCogDist = betterDist, std=globStd[depth])
        tempObjList.append(obj)

tempObjList = sorted(tempObjList, key = lambda x: x.sigmas, reverse = True)
UtilStore(tempObjList, RECLASSIFIED_DIR_LIST(cogDistFuncName))

tempListBySigmas = [x[2] for x in \
    sorted(tempList, key = lambda x: x[0], reverse=True) if x[0] >= 3.0]
tempListByDistDiff = [x[2] for x in \
    sorted(tempList, key = lambda x: x[1], reverse=True) if x[0] >= 3.0]

with open(config.WORK_FILES_DIR() + "Reclassify_" + cogDistFuncName + \
    "_by_sigmas.txt", "w") as f:
    f.write("RECLASSIFICATIONS of genomes\n\n")
    for x in tempListBySigmas:
        f.write(x)

with open(config.WORK_FILES_DIR() + "Reclassify_" + cogDistFuncName + \
    "_by_dist_diff.txt", "w") as f:
    f.write("RECLASSIFICATIONS of genomes\n\n")
    for x in tempListByDistDiff:
        f.write(x)

print("\nReclassified %d good sigmas %d same %d" % (len(tempList),
    tempGoodSigmaCount, tempSameCount))



