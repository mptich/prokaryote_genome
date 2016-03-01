__author__ = 'morel'

import sys
import os
from shared.algorithms.JaccardSuffix import *
from shared.algorithms.string_algo.blosum_matrix import BlosumMatrix
from shared.algorithms.kendall import calculateWeightedKendall
from shared.pyutils.utils import *
import numpy as np

def calculateOrderCorrelation(l):
    """
    :param l - list of pairs
    :return: Kendal Tau correlation
    """
    return calculateWeightedKendall(*(zip(*l)))

def ProcessProteinPair(f, matrix):
    protein1Line = None
    protein2Line = None
    dict12 = {}
    dict21 = {}
    posList12 = []
    posList21 = []

    for i, l in enumerate(f, start = 1):
        l = l.strip()
        if l == "end":
           break
        if i == 1:
            protein1Line = int(l.split('\t')[1])
            length1 = int(l.split('\t')[2])
            continue
        if i == 2:
            protein2Line = int(l.split('\t')[1])
            length2 = int(l.split('\t')[2])
            print("protein1Line %d protein2Line %d" % (protein1Line,
                                                       protein2Line))
            continue

        if protein1Line == protein2Line:
            continue

        if l in ["1", "2"]:
            if (((l=="1") and (length1 < length2)) or
                ((l=="2") and (length1 >= length2))):
                dict = dict12
                posList = posList12
            else:
                dict = dict21
                posList = posList21
            continue

        ll = l.split('\t')
        assert(len(ll) == 4)
        dict[ll[0]] = ll[2]
        posList.append((int(ll[1]), int(ll[3])))

    if (len(dict12) == 0) or (len(dict21) == 0):
        return (None, None, None, None)

    filling12 = len(set(dict12.values())) / float(len(dict12))
    filling21 = len(set(dict21.values())) / float(len(dict21))
    orderCorr12 = calculateOrderCorrelation(posList12)
    orderCorr21 = calculateOrderCorrelation(posList21)
    return (filling12, filling21, orderCorr12, orderCorr21)

    #scoreList = []
    #scoreLossList = []
    #for suff1, suff2 in dict.items():
        #rate, rate1, rate2 = JaccardSuffixDistance(suff1, suff2, matrix)
        #scoreList.append(rate)
        #scoreLossList.append(rate1 - rate)




matrix = BlosumMatrix("shared/algorithms/blosum62.txt")
if len(sys.argv) < 2:
    print ("Missing filename")
    sys.exit(-1)
filename = sys.argv[1]
fillings12 = []
fillings21 = []
orders12 = []
orders21 = []
with open(filename, 'r') as f:
    for i,l in enumerate(f, start=1):
        if l.strip() != "start":
            break
        f12, f21, o12, o21 = ProcessProteinPair(f, matrix)
        if f12 is not None:
            fillings12.append(f12)
            fillings21.append(f21)
            orders12.append(o12)
            orders21.append(o21)
UtilDrawHistogram(fillings12)
UtilDrawHistogram(fillings21)
UtilDrawHistogram(orders12)
UtilDrawHistogram(orders21)
