__author__ = 'morel'

import sys
from shared.algorithms.JaccardSuffix import *
from shared.algorithms.string_algo.blosum_matrix import BlosumMatrix

def ProcessProteinPair(f):
    protein1Line = None
    protein2Line = None
    dict12 = {}
    dict21 = {}

    for i, l in enumerate(f, start = 1):
        l = l.strip()
        if l == "end":
           break
        if i == 1:
            protein1Line = int(l.split('\t')[1])
            continue
        if i == 2:
            protein2Line = int(l.split('\t')[1])
            continue
        if protein1Line == protein2Line:
            continue
        if l == "1":
            dict = dict12
            continue
        if l == "2":
            dict = dict21
            continue
        ll = l.split('\t')
        assert(len(ll) == 2)
        dict[ll[0]] = ll[1]

    if len(dict12) == 0:
        return (None, None None)



matrix = BlosumMatrix("shared/algorithms/blosum62.txt")
if len(sys.argv) < 2:
    print ("Missing filename")
    sys.exit(-1)
filename = sys.argv[1]
scores = []
scoreLosses = []
fills = []
with open(filename, 'r') as f:
    for l in f:
        if l.strip() != "start":
            break
        score, scoreLoss, fill = ProcessProteinPair(f)