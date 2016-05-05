# Place a genome (identified by its directory name) into a Taxonomy
# Classification tree, based on the hierarchical algorithm: average
# COG distances for domain, phylum, etc.

import sys
from taxonomy import *
from collections import defaultdict as DefDict
import common_cogs_method as commonCogsMethod
from shared.algorithms.kendall import calculateWeightedKendall
from shared.pyutils.utils import *
from shared.pyutils.UtilNormDistrib import *
import config
import operator
import math

def prepareForKendall(dirCorrDict, objList, attr):
    dirCorrList = [x for x in dirCorrDict.iteritems() if x[0] in \
        [y.dir for y in objList]]
    corrList = [x[1] for x in sorted(dirCorrList, key = lambda x: x[0])]
    sortedObjList = sorted(objList, key=lambda x: x.dir)
    taxaDistList = [obj.taxaDist for obj in sortedObjList]
    attrList = [getattr(obj, attr) for obj in sortedObjList]
    return (corrList, taxaDistList, attrList)

if len(sys.argv) != 2:
    print "Missing COG distance function name"
    sys.exit(-1)

cogDistFuncName = sys.argv[1]
cogDistFunc = getattr(commonCogsMethod, cogDistFuncName)

cogDict, taxaDict, _, _ = commonCogsMethod.buildCogTaxaDict(3.0)
print ("cogDict len %d, taxaDict len %d" % (len(cogDict), len(taxaDict)))

print("Build taxaType histogram...")
taxaTypeToDirList = DefDict(list)
for dir, taxa in taxaDict.iteritems():
    taxaTypeToDirList[taxa.type].append(dir)
print("taxaTypeToDirList len %d" % len(taxaTypeToDirList))
tempHist = DefDict(int)
for taxaType, l in taxaTypeToDirList.iteritems():
    tempHist[len(l)] += 1
print sorted(tempHist.iteritems(), key=operator.itemgetter(0))

print("Building COG distances...")
cogDist = DefDict(dict)
for ordinal, (dir1, cs1) in enumerate(cogDict.iteritems(), start = 1):
    print("\r%d. %s" % (ordinal, dir1)),
    for dir2, cs2 in cogDict.iteritems():
        cogDist[dir1][dir2] = cogDistFunc(cs1, cs2)

# Build dictionary of taxaType -> list of genome dirs
taxaTypeDirDict = DefDict(list)
print("\nBuilding list of dirs for TaxaTypes...")
for dir, taxa in taxaDict.iteritems():
    taxaTypeDirDict[taxa.type].append(dir)

# Calculate global STD
taxaTypeNormDistDict = DefDict(lambda: None)
print ("Calculating taxaTypeNormDistDict...")
for taxaType, dirList in taxaTypeDirDict.iteritems():
    l = []
    if len(dirList) >= 3:
        for d1 in dirList:
            for d2 in dirList:
                if d1 >= d2:
                    continue
                l.append(cogDist[d1][d2])
        taxaTypeNormDistDict[taxaType] = UtilNormDistrib(list = l)
print("taxaTypeStdDict\n%s" % repr(sorted(taxaTypeNormDistDict.iteritems(),
    key = lambda x: x[1].mean)))
globStd = np.mean([x.std for x in taxaTypeNormDistDict.values()])
globMean = np.mean([x.mean for x in taxaTypeNormDistDict.values()])
print("Global STD %f global mean %f" % (globStd, globMean))

# For each genome, build taxaType -> list of tuples(genome dir,
# COG dist), with the key dir being excluded
dirTaxaTypeDirDict = DefDict(lambda: DefDict(list))
print("\nBuilding list of dirs for TaxaTypes for dir with exclusion...")
for ordinal, dir1 in enumerate(taxaDict.keys(), start = 1):
    print("\r%d. %s" % (ordinal, dir1)),
    for dir2, taxa in taxaDict.iteritems():
        if dir1 == dir2:
            continue
        dirTaxaTypeDirDict[dir1][taxa.type].append((dir2,
            cogDist[dir1][dir2]))

# Build genome reclassifications
tempSameCount = 0
tempImprCount = 0
tempSingleList = []
tempSingleObjList = []
tempList = []
tempObjList = []
print("Build reclassification...")
for dir, taxaTypeDict in dirTaxaTypeDirDict.iteritems():
    oldType = taxaDict[dir].type
    myl = [(np.mean([y[1] for y in x[1]]), x[0]) for x in
        taxaTypeDict.iteritems()]
    bestTaxaType = min(myl)[1]
    if bestTaxaType == oldType:
        tempSameCount += 1
        continue
    l = []
    normDist = UtilNormDistrib(list = [x[1] for x in taxaTypeDict[
        bestTaxaType]])
    newDist = normDist.mean
    taxaDist = bestTaxaType.distance(oldType)
    # Now calculate the mean of distances from dir in its old TaxaType
    if not taxaTypeDict[oldType]:
        distRatio = newDist / globMean
        desc = "%s\n" \
            "Ratio of new distance to average distance: %f\n" \
            "Taxonomy distance between old and new positions: %d\n" \
            "Old classification: %s\n" \
            "New classification: %s\n" \
            "New COG distance: %f\n\n" % \
            (dir, distRatio, taxaDist,
            repr(oldType), repr(bestTaxaType), newDist)
        tempSingleList.append((newDist, desc))
        obj = UtilObject(dir = dir, distRatio = distRatio,
            taxaDist = taxaDist, oldClass = oldType, newClass = bestTaxaType,
            newDist = newDist)
        tempSingleObjList.append(obj)
        continue
    else:
        oldMean = np.mean([x[1] for x in taxaTypeDict[oldType]])
        sigmas = (oldMean - newDist) / globStd
    tempImprCount += 1
    desc = "%s\nReclassification sigmas: %f\n" \
        "Taxonomy distance between old and new positions: %d\n" \
        "Old classification: %s\n" \
        "New classification: %s\n" \
        "Old COG distance: %f\n" \
        "New COG distance: %f\n\n" % \
        (dir, sigmas, taxaDist,
        repr(oldType), repr(bestTaxaType), oldMean, newDist)
    tempList.append((newDist - oldMean, desc))
    obj = UtilObject(dir = dir, sigmas = sigmas,
        taxaDist = taxaDist, oldClass = oldType,
        newClass = bestTaxaType, oldDist = oldMean,
        newDist = newDist)
    tempObjList.append(obj)

tempObjList = sorted(tempObjList, key = lambda x: x.sigmas, reverse = True)
tempSingleObjList = sorted(tempSingleObjList, key = lambda x: x.distRatio)
UtilStore(tempObjList, RECLASSIFIED_DIR_LIST(cogDistFuncName))
UtilStore(tempSingleObjList, RECLASSIFIED_SINGLE_DIR_LIST(cogDistFuncName))

tempList = [x[1] for x in sorted(tempList)]
tempSingleList = [x[1] for x in sorted(tempSingleList)]

with open(config.WORK_FILES_DIR() + "Reclassify_" + cogDistFuncName + \
    ".txt", "w") as f:
    f.write("RECLASSIFICATIONS of "\
        "genomes with neighbors in the same family\n\n")
    for x in tempList:
        f.write(x)
    f.write("\n\n\n")
    f.write("RECLASSIFICATIONS of genomes with their own taxonomy family\n")
    f.write("All of them are moved into another, " \
        "already populated family\n")
    f.write("Likelihood of the reclassification is higher for genomes with\n")
    f.write("low ratio of new distance to average distance\n\n")
    for x in tempSingleList:
        f.write(x)

print("\nReclassified %d single %d same %d" % (tempImprCount,
    len(tempSingleList), tempSameCount))

dirCorrDict = UtilLoad(GENOME_CORR_DICT())
dirCorrSelList, taxaDistList, attrList = prepareForKendall(dirCorrDict,
    tempObjList, "sigmas")
correlation = calculateWeightedKendall(taxaDistList, dirCorrSelList)
attrCorrelation = calculateWeightedKendall(attrList, dirCorrSelList)
dirCorrSelList, taxaDistList, attrList = prepareForKendall(dirCorrDict,
    tempSingleObjList, "distRatio")
correlationSingle = calculateWeightedKendall(taxaDistList, dirCorrSelList)
attrCorrelationSingle = calculateWeightedKendall(attrList, dirCorrSelList)

print("Correlation for reclass %f, for singles %f" % (correlation,
    correlationSingle))
print("Attr correlation for reclass %f, for singles %f" % (attrCorrelation,
    attrCorrelationSingle))


