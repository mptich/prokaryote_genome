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
from shared.pyutils.distance_matrix import *

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

print("reading COG instance set...")
with open(COG_INST_SET(), 'r') as fset:
    cogSet = json.load(fset, object_hook = UtilJSONDecoderDictToObj)
print("Read %d COG instances" % len(cogSet))

print ("processing COG instance set...")
cogDict = {}
cogCounter = 0 # converting name -> number for speed
cogNameToIdDict = {}
for cogInst in cogSet:
    dir = cogInst.getDir()
    cs = cogDict.get(dir, set())
    name = cogInst.name
    if name not in cogNameToIdDict:
        cogNameToIdDict[name] = cogCounter
        cogCounter += 1
    cs.add(cogNameToIdDict[name])
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

print("Building COG distances...")
cogDist = {}
cogWeights = {}
cogWeightsReverse = {}
for dir1, cs1 in cogDict.items():
    if dir1 not in validDirSet:
        continue
    dictDirCogDist = {}
    cogDist[dir1] = dictDirCogDist
    dictDirWeights = {}
    cogWeights[dir1] = dictDirWeights
    dictDirWeightsReverse = {}
    cogWeightsReverse[dir1] = dictDirWeightsReverse
    for dir2, cs2 in cogDict.items():
        if dir2 not in validDirSet:
            continue
        dictDirCogDist[dir2] = commonCogsDist(cs1, cs2)
        dictDirWeights[dir2] = commonCogsWeight(cs1, cs2)
        dictDirWeightsReverse[dir2] = commonCogsWeightReverse(cs1, cs2)

print("Building Taxonomy distances...")
taxDist = {}
for dir1, taxa1 in taxaDict.items():
    if dir1 not in validDirSet:
        continue
    dict = {}
    taxDist[dir1] = dict
    for dir2, taxa2 in taxaDict.items():
        if dir2 not in validDirSet:
            continue
        d = taxa1.distance(taxa2)
        d += ((random.randrange(10000) - 5000) / 15000.)
        dict[dir2] = d

cogDistMat = DistanceMatrix(doubleDict=cogDist)
print("Got cogDistMat")
cogWeightsMat = DistanceMatrix(doubleDict=cogWeights)
print("Got cogWeightsMat")
cogWeightsReverseMat = DistanceMatrix(doubleDict=cogWeightsReverse)
print("Got cogWeightsReverseMat")
taxDistMat = DistanceMatrix(doubleDict=taxDist)
print("Got taxDistMat")

mean, std, corrList = distanceMatrixCorrelation(taxDistMat, cogDistMat,
                                         cogWeightsMat)
print("Weighted correlation: mean %f std %f" % (mean, std))
print "Worst correlations: ", corrList[:10]
mean, std, corrList = distanceMatrixCorrelation(taxDistMat, cogDistMat,
                                         cogWeightsReverseMat)
print("Reverse weighted correlation: mean %f std %f" % (mean, std))
print "Worst correlations: ", corrList[:10]
mean, std, corrList = distanceMatrixCorrelation(taxDistMat, cogDistMat, None)
print("Unweighted correlation: mean %f std %f" % (mean, std))
print "Worst correlations: ", corrList[:10]






