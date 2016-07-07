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

def cogSetWeight(cogSet, cogFreq, cogReg):
    return sum([1. / (cogFreq[cogName] + cogReg) for cogName in cogSet])


def commonCogsDist(dir1, dir2, cogDict):
    cs1 = cogDict[dir1]
    cs2 = cogDict[dir2]
    commonSet = cs1 & cs2
    l = len(commonSet) + 1
    return math.log(float(len(cs1) + 1) * (len(cs2) + 1) / (l * l))

def commonCogsDistAdj(dir1, dir2, cogDict, cogFreq):
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
    cs1 = cogDict[dir1]
    cs2 = cogDict[dir2]
    commonSet = cs1 & cs2
    c = float(len(commonSet))
    na = float(len(cs1))
    nb = float(len(cs2))
    t = float(len(cogFreq))
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

MaxRegDist = 10
ExpMaxRegDist = math.exp(MaxRegDist)
def commonCogsDistReg(dir1, dir2, cogDict, cogWeightDict, reg):
    """
    Calculates regularized distances based on common COGs
    :param cs1: 1st set of cog names
    :param cs2: 2nd set of cog names
    :return: distance, regularized
    """
    cs1 = cogDict[dir1]
    cs2 = cogDict[dir2]
    commonSetWeight = cogWeightDict[dir1][dir2]
    lenList = [len(x) for x in [cs1, cs2]]
    weightList = [cogWeightDict[dir1][dir1], cogWeightDict[dir2][dir2]]
    minSetWeight = weightList[lenList.index(min(lenList))] + reg
    return math.log((ExpMaxRegDist + 1.0) * minSetWeight /
        (ExpMaxRegDist * commonSetWeight + minSetWeight))


def buildCogTaxaDict(showCogFreqHist = False):

    print("Reading cogDict...")
    cogDict = UtilLoad(COG_DICT())

    print("Building COG frequncies...")
    cogFreq = DefDict(int)
    for dir, cogs in cogDict.iteritems():
        for cname in cogs:
            cogFreq[cname] += 1

    if showCogFreqHist:
        print("Sowing cogFreq histogram...")
        UtilDrawHistogram(cogFreq.values(), show = True)
        print("Showing genome specificity histogram...")
        UtilDrawHistogram([x for y in cogWeightDict.values() \
            for x in y.values()], show = True)

    print("reading taxa dictionary...")
    taxaDict = UtilLoad(PROK_TAXA_DICT())
    print("Read %d organisms" % len(taxaDict))

    temp = taxaDict.keys()
    for dir in temp:
        if dir not in cogDict:
            del taxaDict[dir]
    temp = cogDict.keys()
    for dir in temp:
        if dir not in taxaDict:
            del cogDict[dir]
    print("Valid set contains %d organisms" % len(cogDict))

    return (cogDict, cogFreq, taxaDict)


def buildCogDistances(cogDict, cogFreq, cogReg, genReg):

    print("Building cogWeightsDict...")
    cogWeightDict = DefDict(dict)
    for ind, (dir1, cogs1) in enumerate(cogDict.iteritems(), start=1):
        print("\r%d. %s" % (ind, dir1)),
        for dir2, cogs2 in cogDict.iteritems():
            cogWeightDict[dir1][dir2] =\
                cogSetWeight(cogs1 & cogs2, cogFreq, cogReg)

    print("\nBuilding COG distances...")
    cogDist = DefDict(dict)
    for ordinal, dir1 in enumerate(cogDict, start = 1):
        print("\r%d. %s" % (ordinal, dir1)),
        for dir2 in cogDict:
            cogDist[dir1][dir2] = commonCogsDistReg(dir1, dir2, cogDict,
                cogWeightDict, genReg)

    print("\nStoring COG distance dictionary...")
    UtilStore(cogDist, COG_DIST_DICT(commonCogsDistReg.__name__))


if __name__ == "__main__":

    cogDict, cogFreq, taxaDict = buildCogTaxaDict()
    buildCogDistances(cogDict, cogFreq, 50., 0.5)

    print("Reading COG distances...")
    cogDist = UtilLoad(COG_DIST_DICT(commonCogsDistReg.__name__))

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
    taxDistMat = DistanceMatrix(doubleDict=taxDist)
    print("Got taxDistMat")
    # randDistMat = DistanceMatrix(doubleDict=randDist)
    # print("Got randDistMat")

    mean, std, corrList, comp = distanceMatrixCorrelation(taxDistMat,
        cogDistMat, None, True)
    print("Correlation: mean %f std %f" % (mean, std))
    print "Components ", comp
    print "Worst correlations: ", corrList[:10]
    print "Best correlations: ", corrList[-10:], "\n"
    #UtilDrawHistogram(inputList = [x[1] for x in corrList], show = False)

    UtilStore(dict(corrList), GENOME_CORR_DICT())

    # mean, std, corrList, _ = distanceMatrixCorrelation(taxDistMat,
        # randDistMat, None, False)
    # print("Unweighted totally random correlation: mean %f std %f" %
        # (mean, std))







