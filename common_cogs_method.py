# This module analyzes correlation between distances between genomes
# defined as share of common COGs, and distances between genomes deifned by
# their taxonomy.
# If genome A has Na COGs, and genome B got Nb cogs, and they have N COGs in
# common, then their distance in terms of COGs is defined as (Na * Nb) /
# ((N+1) * (N+1))
# Distance between COGs in terms of their taxonomy is defined in
# Taxa.distance() method (file taxonomy.py)

from taxonomy import *
import random
import math
from collections import defaultdict as DefDict
from shared.pyutils.distance_matrix import *

# Set all COG names
cogNameSet = set()

def commonCogsDist(cs1, cs2):
    """
    Calculates distances based on common COGs
    :param cs1: 1st set of cog names
    :param cs2: 2nd set of cog names
    :return: distance, based on the equation (len(cs1) * len(cs2)) / ((N+1)
    * (N+1)), where N is the number of common COG names
    """
    commonSet = cs1 & cs2
    l = len(commonSet) + 1
    return float(len(cs1)) * len(cs2) / (l * l)

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

def commonCogsDistNonRand(cs1, cs2):
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
            d = (Na*Nb) / ((EC+1)^2)
    """
    commonSet = cs1 & cs2
    c = float(len(commonSet))
    na = float(len(cs1))
    nb = float(len(cs2))
    t = float(len(cogNameSet))
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
    return (na * nb) / (ec * ec)


print("reading COG instance set...")
with open(COG_INST_SET(), 'r') as fset:
    cogSet = json.load(fset, object_hook = UtilJSONDecoderDictToObj)
print("Read %d COG instances" % len(cogSet))

print ("processing COG instance set...")
cogDict = {}
for cogInst in cogSet:
    dir = cogInst.getDir()
    cs = cogDict.get(dir, set())
    name = cogInst.getName()
    cs.add(name)
    cogDict[dir] = cs
print("Got %d organisms with COGS" % len(cogDict))

print("reading taxa dictionary...")
with open(PROK_TAXA_DICT(), 'r') as fdict:
    taxaDict = json.load(fdict, object_hook = UtilJSONDecoderDictToObj)
print("Read %d organisms" % len(taxaDict))

validDirSet = set()
for dir in taxaDict:
    if dir in cogDict:
        validDirSet.add(dir)
print("Valid set contains %d organisms" % len(validDirSet))

# Build COG name set so it includes only COG from genomes in validDirSet
for dir in validDirSet:
    cogNameSet |= cogDict[dir]
print("Got %d total COGs in valid genomes" % len(cogNameSet))

# Building random distance for comparison
print("Building random distances...")
randDist = DefDict(dict)
for dir1 in validDirSet:
    for dir2 in validDirSet:
        randDist[dir1][dir2] = random.random()

print("Building COG distances...")
cogDist = DefDict(dict)
cogDistNonRand = DefDict(dict)
cogWeights = DefDict(dict)
cogWeightsReverse = DefDict(dict)
for ordinal, (dir1, cs1) in enumerate(cogDict.items(), start = 1):
    if dir1 not in validDirSet:
        continue
    print("\r%d. %s" % (ordinal, dir1)),
    for dir2, cs2 in cogDict.items():
        if dir2 not in validDirSet:
            continue
        cogDist[dir1][dir2] = commonCogsDist(cs1, cs2)
        cogDistNonRand[dir1][dir2] = commonCogsDistNonRand(cs1, cs2)
        cogWeights[dir1][dir2] = commonCogsWeight(cs1, cs2)
        cogWeightsReverse[dir1][dir2] = commonCogsWeightReverse(cs1, cs2)

print("\nBuilding Taxonomy distances...")
taxDist = DefDict(dict)
for dir1, taxa1 in taxaDict.items():
    if dir1 not in validDirSet:
        continue
    for dir2, taxa2 in taxaDict.items():
        if dir2 not in validDirSet:
            continue
        d = taxa1.distance(taxa2)
        taxDist[dir1][dir2] = d

cogDistMat = DistanceMatrix(doubleDict=cogDist)
print("Got cogDistMat")
cogDistNonRandMat = DistanceMatrix(doubleDict=cogDistNonRand)
print("Got cogDistNonRandMat")
cogWeightsMat = DistanceMatrix(doubleDict=cogWeights)
print("Got cogWeightsMat")
cogWeightsReverseMat = DistanceMatrix(doubleDict=cogWeightsReverse)
print("Got cogWeightsReverseMat")
taxDistMat = DistanceMatrix(doubleDict=taxDist)
print("Got taxDistMat")
randDistMat = DistanceMatrix(doubleDict=randDist)
print("Got randDistMat")

mean, std, corrList, comp = distanceMatrixCorrelation(taxDistMat, cogDistMat,
                                         cogWeightsMat, True)
print("Weighted correlation: mean %f std %f" % (mean, std))
print "Components ", comp
print "Worst correlations: ", corrList[:10], "\n"

mean, std, corrList, comp = distanceMatrixCorrelation(taxDistMat, cogDistMat,
                                         cogWeightsReverseMat, True)
print("Reverse weighted correlation: mean %f std %f" % (mean, std))
print "Components ", comp
print "Worst correlations: ", corrList[:10], "\n"

mean, std, corrList, comp = distanceMatrixCorrelation(taxDistMat, cogDistMat,
                                                None, True)
print("Unweighted correlation: mean %f std %f" % (mean, std))
print "Components ", comp
print "Worst correlations: ", corrList[:10], "\n"

mean, std, corrList, comp = distanceMatrixCorrelation(taxDistMat,
                                                  cogDistNonRandMat,
                                                None, True)
print("Unweighted non random correlation: mean %f std %f" % (mean, std))
print "Components ", comp
print "Worst correlations: ", corrList[:10], "\n"

mean, std, corrList, _ = distanceMatrixCorrelation(taxDistMat, randDistMat,
                                                None, False)
print("Unweighted totally random correlation: mean %f std %f" % (mean, std))






