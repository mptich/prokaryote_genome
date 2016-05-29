# This module analyzes correlation between distances between genomes
# defined as share of common COGs, and distances between genomes deifned by
# their taxonomy.
# If genome A has Na COGs, and genome B got Nb cogs, and they have N COGs in
# common, then their distance in terms of COGs is defined as (Na * Nb) /
# ((N+1) * (N+1))
# Distance between COGs in terms of their taxonomy is defined in
# TaxaType.distance() method (file taxonomy.py)

from taxonomy import *
import random
import math
import sys
import operator
from collections import defaultdict as DefDict
from shared.pyutils.utils import *
from shared.pyutils.distance_matrix import *
from shared.algorithms.kendall import calculateWeightedKendall

# COG name -> list of its lengthes
cogLengthDict = DefDict(list)

# Total summary statistical weight of all COGs
totalCogWeight = 0.

def commonCogsDist(cs1, cs2):
    """
    Calculates distances based on common COGs
    :param cs1: 1st set of cog names
    :param cs2: 2nd set of cog names
    :return: distance, based on the equation ln((len(cs1) + 1) * (len(cs2) + 1)
        / ((N+1) * (N+1))), where N is the number of common COG names
    """
    commonSet = cs1 & cs2
    l = len(commonSet) + 1
    return math.log(float(len(cs1) + 1) * (len(cs2) + 1) / (l * l))

def commonCogsWeight(cs1, cs2):
    """
    Calculates weight of the distance based on common COGs
    :param cs1: 1st set of cog names
    :param cs2: 2nd set of cog names
    :return: weight, based on the equation sqrt(len(cs1) * len(cs2))
    """
    return math.sqrt(float(len(cs1)) * len(cs2))

def commonCogsWeightReverse(cs1, cs2):
    """
    Calculates weight of the distance based on common COGs
    :param cs1: 1st set of cog names
    :param cs2: 2nd set of cog names
    :return: weight, based on the equation 1 / sqrt((len(cs1)+1) * (len(
            cs2)+1))
    """
    return 1. / math.sqrt(float(len(cs1) + 1) * (len(cs2) + 1))

def commonCogsDistAdj(cs1, cs2):
    """
    Calculates weight of the distance based on common COGs
    :param cs1: 1st set of cog names
    :param cs2: 2nd set of cog names
    :return: Weight takes into account only non-random common
            COGs count. Math behind this calculation is described in the
            following article:
            TODO
            Equation is as follows:
            if Na - # of COGs in genome A (len(cs1)),
            Nb - number of COGs in genome B (len(cs2)),
            C - measured number of common COGs,
            T - total number of considered COGs,
            then
            EC = (T*C - Na*Nb) / (T+C-Na-Nb) - expected # of common COGs (
            take it as 0 if EC < 0), and the distance between genomes will be
            d = ln((Na+1)*(Nb+1) / ((EC+1)^2))
    """
    commonSet = cs1 & cs2
    c = float(len(commonSet))
    na = float(len(cs1))
    nb = float(len(cs2))
    t = float(len(cogLengthDict))
    assert(t >= na)
    assert(t >= nb)
    denom = t + c - na - nb
    if denom:
        ec = (t * c - na * nb) / denom
    else:
        ec = -1.
    if ec < 0.:
        ec = 0.
    ec += 1.
    return math.log((na + 1) * (nb + 1) / (ec * ec))


def buildCogTaxaDict(cogLengthFilter, standAloneList = None):
    if standAloneList is None:
        standAloneList = []

    print("reading COG instance list...")
    cogList = UtilLoad(COG_INST_LIST())
    print("Read %d COG instances" % len(cogList))

    print ("Building cogLengthDict...")
    for cogInst in cogList:
        cogLengthDict[cogInst.name].append(cogInst.len)
    print("COGs read from file: %d" % len(cogLengthDict))

    print ("Translating cogLengthDict...")
    temp = cogLengthDict.keys()
    for cogName in temp:
        lengthList = cogLengthDict[cogName]
        cogLengthDict[cogName] = (np.mean(lengthList), np.std(lengthList))

    print ("Building cogDict...")
    cogDict = DefDict(set)
    validCogInstances = 0
    for cogInst in cogList:
        name = cogInst.name
        assert(name in cogLengthDict)
        length = cogInst.len
        lenMean, lenStd = cogLengthDict[name]
        if (length < (lenMean - cogLengthFilter * lenStd)) or \
                (length > (lenMean + cogLengthFilter * lenStd)):
            continue
        validCogInstances += 1
        dir = cogInst.dir
        cogDict[dir].add(name)
    print("Got %d organisms with COGS" % len(cogDict))
    print("Read %d COG instances, selected %d out of them" %
        (len(cogList), validCogInstances))

    print("reading taxa dictionary...")
    taxaDict = UtilLoad(PROK_TAXA_DICT())
    print("Read %d organisms" % len(taxaDict))

    standAloneCogDict = {}
    for dir in standAloneList:
        if dir not in cogDict:
            raise ValueError("Bad genome dir %s" % dir)
        standAloneCogDict[dir] = cogDict[dir]

    standAloneTaxaDict = {}
    for dir in standAloneList:
        if dir not in taxaDict:
            standAloneTaxaDict[dir] = "Unknown"
        else:
            standAloneTaxaDict[dir] = taxaDict[dir]

    temp = taxaDict.keys()
    for dir in temp:
        if dir not in cogDict:
            del taxaDict[dir]
    temp = cogDict.keys()
    for dir in temp:
        if dir not in taxaDict:
            del cogDict[dir]
    print("Valid set contains %d organisms" % len(cogDict))

    return (cogDict, taxaDict, standAloneCogDict, standAloneTaxaDict)


if __name__ == "__main__":

    if len(sys.argv) == 1:
        cogLengthFilter = 1.0
    else:
        cogLengthFilter = float(sys.argv[1])
    print ("Using cogLengthFilter %f" % cogLengthFilter)

    cogDict, taxaDict, _, _ = buildCogTaxaDict(cogLengthFilter)

    # Building random distance for comparison
    # print("Building random distances...")
    # randDist = DefDict(dict)
    # for dir1 in cogDict:
        # for dir2 in cogDict:
            # randDist[dir1][dir2] = random.random()

    print("Building COG distances...")
    cogDist = DefDict(dict)
    cogDistStat = DefDict(dict)
    for ordinal, (dir1, cs1) in enumerate(cogDict.items(), start = 1):
        print("\r%d. %s" % (ordinal, dir1)),
        for dir2, cs2 in cogDict.items():
            cogDist[dir1][dir2] = commonCogsDist(cs1, cs2)
            cogDistStat[dir1][dir2] = commonCogsDistAdj(cs1, cs2)

    print("\nBuilding Taxonomy distances...")
    taxDist = DefDict(dict)
    for dir1, taxa1 in taxaDict.items():
        for dir2, taxa2 in taxaDict.items():
            d = taxa1.distance(taxa2)
            taxDist[dir1][dir2] = d

    print("Building dict of taxonomy dist counts...")
    genTaxDistCntDict = DefDict(lambda: [0] * (TaxaType.maxDistance() + 1))
    for dir, tdd in taxDist.items():
        for d in tdd.values():
            genTaxDistCntDict[dir][d] += 1
    UtilStore(genTaxDistCntDict, GENOME_TAX_DIST_CNT_DICT())
    ttTaxDistCntDict = {}
    for dir, l in genTaxDistCntDict.items():
        ttTaxDistCntDict[taxaDict[dir].type.key] = l
    UtilStore(ttTaxDistCntDict, TAXTYPE_TAX_DIST_CNT_DICT())

    cogDistMat = DistanceMatrix(doubleDict=cogDist)
    print("Got cogDistMat")
    cogDistStatMat = DistanceMatrix(doubleDict=cogDistStat)
    print("Got cogDistNonRandMat")
    taxDistMat = DistanceMatrix(doubleDict=taxDist)
    print("Got taxDistMat")
    # randDistMat = DistanceMatrix(doubleDict=randDist)
    # print("Got randDistMat")

    mean, std, corrList, comp = distanceMatrixCorrelation(taxDistMat,
        cogDistMat, None, True)
    print("Unweighted correlation: mean %f std %f" % (mean, std))
    print "Components ", comp
    print "Worst correlations: ", corrList[:10]
    print "Best correlations: ", corrList[-10:], "\n"
    UtilDrawHistogram(inputList = [x[1] for x in corrList], show = False)

    mean, std, corrList, comp = distanceMatrixCorrelation(taxDistMat,
                                                      cogDistStatMat,
                                                    None, True)
    print("Unweighted non random correlation: mean %f std %f" % (mean, std))
    print "Components ", comp
    print "Worst correlations: ", corrList[:10]
    print "Best correlations: ", corrList[-10:], "\n"
    UtilDrawHistogram(inputList = [x[1] for x in corrList], show = True)

    UtilStore(dict(corrList), GENOME_CORR_DICT())

    # mean, std, corrList, _ = distanceMatrixCorrelation(taxDistMat,
        # randDistMat, None, False)
    # print("Unweighted totally random correlation: mean %f std %f" %
        # (mean, std))







