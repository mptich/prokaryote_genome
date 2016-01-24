# This module analyzes correlation between distances between genomes
# defined as share of common COGs, and distances between genomes deifned by
# their taxonomy.
# If genome A has Na COGs, and genome B got Nb cogs, and they have N COGs in
# common, then their distance in terms of COGs is defined as (Na * Nb) /
# ((N+1) * (N+1))
# Distance between COGs in terms of their taxonomy is defined in
# Taxa.distance() method (file taxonomy.py)

from import_proxy import *
from filedefs import *
import numpy as np
from taxonomy import *
import random
import sys

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
for dir1, cs1 in cogDict.items():
    if dir1 not in validDirSet:
        continue
    dict = {}
    cogDist[dir1] = dict
    for dir2, cs2 in cogDict.items():
        if dir2 not in validDirSet:
            continue
        dict[dir2] = commonCogsDist(cs1, cs2)

print("Building Taxonomy distances...")
histTaxDist = [0] * (Taxa.maxDistance() + 1)

cogDistLists = []
for i in range(Taxa.maxDistance() + 1):
    cogDistLists.append([])
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
        histTaxDist[d] += 1
        cogDistLists[d].append(cogDist[dir1][dir2])
        d += ((random.randrange(10000) - 5000) / 15000.)
        dict[dir2] = d
print("Histogram of taxonomy distances: %s" %
      ', '.join(map(str, histTaxDist)))
print("COG distannce distribution %s" % ', '.join(map(str, [(np.mean(x),
        np.std(x)) for x in cogDistLists])))

print("Comparing taxonomy and common COG distances...")
distList = []
weights = [1.] * len(validDirSet)
for dir in validDirSet:
    taxDict = taxDist[dir]
    cogDict = cogDist[dir]
    taxList = [taxDict[x] for x in validDirSet]
    #weights = [1/(x+1.) for x in taxList]
    cogList = [cogDict[x] for x in validDirSet]
    d = calculateWeightedKendall(cogList, taxList, weights)
    distList.append(d)

print("Dist correlation: mean %f std %f" % (np.mean(distList), np.std(
    distList)))






