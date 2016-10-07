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
import os.path
import operator
from collections import defaultdict as DefDict
from shared.pyutils.utils import *
from shared.pyutils.distance_matrix import *
from shared.algorithms.kendall import calculateWeightedKendall
from scipy.optimize import anneal

CogDistOptimalParams = \
    {"cogReg" : 5.44122751, "genReg" : -5.85405896, "mixReg" : 0.17919745}

# Unit of quantization of COG weight regularization
COG_REG_STEP = 0.5
COG_REG_LOWER = 0.
COG_REG_STEP_COUNT = 16
CogRegExpSteps = [math.exp(COG_REG_LOWER + float(i) * COG_REG_STEP) for i \
    in range(0, COG_REG_STEP_COUNT+1)]

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
def commonCogsDistReg(dir1, dir2, cogDict, cogWeightDict, expGenReg, mixReg):
    """
    Calculates regularized distances based on common COGs
    :param cs1: 1st set of cog names
    :param cs2: 2nd set of cog names
    :return: distance, regularized
    """
    cs1 = cogDict[dir1]
    cs2 = cogDict[dir2]
    commonSetWeight = cogWeightDict[dir1][dir2] * (1. + mixReg)
    lenList = [len(x) for x in [cs1, cs2]]
    minLenIndex = lenList.index(min(lenList))
    maxLenIndex = 1 - minLenIndex
    weightList = [cogWeightDict[dir1][dir1], cogWeightDict[dir2][dir2]]
    minSetWeight = weightList[minLenIndex] + expGenReg
    maxSetWeight = weightList[maxLenIndex] + expGenReg
    setWeight = minSetWeight + mixReg * maxSetWeight
    return math.log((ExpMaxRegDist + 1.0) * setWeight /
        (ExpMaxRegDist * commonSetWeight + setWeight))


def buildCogTaxaDict(noWeights = False, showCogFreqHist = False,
    interpolationRange = None):

    print("reading taxa dictionary...")
    taxaDict = UtilLoad(PROK_TAXA_DICT())
    print("Read %d organisms" % len(taxaDict))

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

    temp = taxaDict.keys()
    for dir in temp:
        if dir not in cogDict:
            del taxaDict[dir]
    temp = cogDict.keys()
    for dir in temp:
        if dir not in taxaDict:
            del cogDict[dir]
    print("Valid set contains %d organisms" % len(cogDict))

    print("\nBuilding Taxonomy distances...")
    taxDist = DefDict(dict)
    for dir1, taxa1 in taxaDict.items():
        for dir2, taxa2 in taxaDict.items():
            d = taxa1.distance(taxa2)
            taxDist[dir1][dir2] = d

    # Optimization
    if noWeights:
        return (cogDict, None, taxaDict, taxDist)

    fname = COG_WEIGHTS_DICT_LIST()
    if os.path.isfile(fname):
        print("Loading cogWeightDictList...")
        cogWeightDictList = UtilLoad(fname, progrIndPeriod=100)
    else:
        print("Building cogWeightsDict...")
        cogWeightDictList = [DefDict(dict) for i \
            in range(0, COG_REG_STEP_COUNT+1)]
        if not interpolationRange:
            interpolationRange = range(0, COG_REG_STEP_COUNT+1)
        for i in interpolationRange:
            expCogReg = math.exp(COG_REG_LOWER + float(i) * COG_REG_STEP)
            print("\nexpCogReg %f" % expCogReg)
            cogWeightDict = cogWeightDictList[i]
            for ind, (dir1, cogs1) in enumerate(cogDict.iteritems(), start=1):
                print("\r%d.%d. %s" % (i, ind, dir1)),
                for dir2, cogs2 in cogDict.iteritems():
                    cogWeightDict[dir1][dir2] = \
                        cogSetWeight(cogs1 & cogs2, cogFreq, expCogReg)
            print
        UtilStore(cogWeightDictList, fname)

    return (cogDict, cogWeightDictList, taxaDict, taxDist)


def calculateCogRegInt(cogReg):
    cogRegInt = int((cogReg - COG_REG_LOWER) / COG_REG_STEP)
    assert(cogRegInt <= COG_REG_STEP_COUNT)
    if cogRegInt == COG_REG_STEP_COUNT:
        cogRegInt = COG_REG_STEP_COUNT-1
    return cogRegInt


def buildCogDistances(cogDict, cogWeightDictList, cogReg, genReg, mixReg):

    expGenReg = math.exp(genReg)

    cogRegInt = calculateCogRegInt(cogReg)

    expCogReg = math.exp(cogReg)
    fraction = (expCogReg - CogRegExpSteps[cogRegInt]) / \
        (CogRegExpSteps[cogRegInt+1] - CogRegExpSteps[cogRegInt])
    print("Building COG weights interpolation from %d fraction %f" %
          (cogRegInt, fraction))
    cogWeightDictLow = cogWeightDictList[cogRegInt]
    cogWeightDictUpper = cogWeightDictList[cogRegInt+1]
    cogWeightDict = DefDict(dict)
    for dir1, dd in cogWeightDictLow.iteritems():
        for dir2, wl in dd.iteritems():
            wu = cogWeightDictUpper[dir1][dir2]
            cogWeightDict[dir1][dir2] = wl + (wu - wl) * fraction

    print("\nBuilding COG distances...")
    cogDist = DefDict(dict)
    for ordinal, dir1 in enumerate(cogDict, start = 1):
        print("\r%d. %s" % (ordinal, dir1)),
        for dir2 in cogDict:
            cogDist[dir1][dir2] = commonCogsDistReg(dir1, dir2, cogDict,
                cogWeightDict, expGenReg, mixReg)

    return cogDist

def calculateCorrelation(cogDist, taxDist):
    corrList = []
    for dir, cogDirDist in cogDist.iteritems():
        taxDirDist = taxDist[dir]
        corrList.append(calculateWeightedKendall( \
            [taxDirDist[x] for x in cogDirDist.keys()], cogDirDist.values()))

    mean = np.mean(corrList)
    std = np.std(corrList, ddof = 1.)
    print("Result: mean %f std %f" % (mean, std))
    return(mean, std)

bestCorr = 0.909499
bestParamVector = [ 5.44122751, -5.85405896,  0.17919745]
def optimizingFunction(taxDist, cogDict, cogWeightDictList, lowerBounds,
    upperBounds, paramVector):
    global bestCorr, bestParamVector
    cogReg, genReg, mixReg = paramVector
    print("Optimizing function cogReg=%f genReg=%f mixReg=%f" %
        (cogReg, genReg, mixReg))
    for ind, val in enumerate(paramVector):
        if (val < lowerBounds[ind]) or (val > upperBounds[ind]):
            print("Out of bounds index %d" % ind)
            return 1000.
    cogDist = buildCogDistances(cogDict, cogWeightDictList,
        cogReg, genReg, mixReg)
    corr, std = calculateCorrelation(cogDist, taxDist)
    print("CORRELATION: %f STD: %f" % (corr, std))
    if corr > bestCorr:
        bestCorr = corr
        bestParamVector = paramVector
    print("BEST SO FAR: corr %f paramVector %s" % (bestCorr,
        repr(bestParamVector)))
    return 1. / (corr - 1.)


def findOptimum(cogDict, cogWeightDictList, taxDist):
    """
    Finds values of cogReg, genReg, mixReg achieving maximum
    correlation between cogDist and taxDist
    :return: None
    """

    lowerBounds = [COG_REG_LOWER, -7., 0.12]
    upperBounds = [COG_REG_LOWER + COG_REG_STEP * COG_REG_STEP_COUNT, -1.0,
        0.26]
    cycleCount = 0

    while True:
        initParamVector = bestParamVector
        for i in range(len(lowerBounds)):
            lowerBounds[i] += (bestParamVector[i] - lowerBounds[i]) * 0.2
            upperBounds[i] += (bestParamVector[i] - upperBounds[i]) * 0.2

        cycleCount += 1
        print "CYCLE ", cycleCount
        print "initParamVector ", initParamVector
        print "BOUNDS ", lowerBounds, upperBounds
        func = UtilCaller(optimizingFunction, taxDist, cogDict,
            cogWeightDictList, lowerBounds, upperBounds)
        result = anneal(func, initParamVector, schedule='boltzmann',
            full_output=True, maxiter=80, lower=lowerBounds,
            upper=upperBounds, disp=True)
        print result


if __name__ == "__main__":

    """
    Takes the following command line options:
    buildWeights - building COG weights dictionary
    optimize - find optimal parameters for COG weights
    store <cogReg> <genReg> <mixReg> - stores COG distance dictionary
    distCounts - buils taxonomy distance dictionaries
    """

    if (len(sys.argv) == 2) and (sys.argv[1] == "buildWeights"):
        cogDict, cogWeightDictList, taxaDict, taxDist = buildCogTaxaDict()
        # Done
        sys.exit(0)

    if (len(sys.argv) == 2) and (sys.argv[1] == "optimize"):
        cogDict, cogWeightDictList, taxaDict, taxDist = buildCogTaxaDict()
        findOptimum(cogDict, cogWeightDictList, taxDist)
        sys.exit(0)

    if (len(sys.argv) == 5) and (sys.argv[1] == "store"):
        cogReg = float(sys.argv[2])
        cogRegInt = calculateCogRegInt(cogReg)
        cogDict, cogWeightDictList, taxaDict, taxDist = \
            buildCogTaxaDict(interpolationRange=range(cogRegInt, cogRegInt+2))
        cogDist = buildCogDistances(cogDict, cogWeightDictList,
            cogReg, float(sys.argv[3]), float(sys.argv[4]))
        corr, std = calculateCorrelation(cogDist, taxDist)
        print("CORRELATION: %f STD: %f" % (corr, std))
        print("\nStoring COG distance dictionary...")
        UtilStore(cogDist, COG_DIST_DICT())
        sys.exit(0)

    if (len(sys.argv) == 2) and (sys.argv[1] == "optimalStore"):
        cogReg = CogDistOptimalParams["cogReg"]
        cogRegInt = calculateCogRegInt(cogReg)
        cogDict, cogWeightDictList, taxaDict, taxDist = \
            buildCogTaxaDict(interpolationRange=range(cogRegInt, cogRegInt+2))
        cogDist = buildCogDistances(cogDict, cogWeightDictList,
            **CogDistOptimalParams)
        corr, std = calculateCorrelation(cogDist, taxDist)
        print("CORRELATION: %f STD: %f" % (corr, std))
        print("\nStoring COG distance dictionary...")
        UtilStore(cogDist, COG_DIST_DICT())
        sys.exit(0)

    if (len(sys.argv) == 2) and (sys.argv[1] == "distCounts"):
        print("Building dict of taxonomy dist counts...")
        _, _, taxaDict, taxDist = \
            buildCogTaxaDict(noWeights = True)
        genTaxDistCntDict = DefDict(lambda: [0] *
            (TaxaType.maxDistance() + 1))
        for dir, tdd in taxDist.items():
            for d in tdd.values():
                genTaxDistCntDict[dir][d] += 1
        UtilStore(genTaxDistCntDict, GENOME_TAX_DIST_CNT_DICT())
        ttTaxDistCntDict = {}
        for dir, l in genTaxDistCntDict.items():
            ttTaxDistCntDict[taxaDict[dir].type.key] = l
        UtilStore(ttTaxDistCntDict, TAXTYPE_TAX_DIST_CNT_DICT())
        sys.exit(0)

    print("WRONG COMMAND LINE")







