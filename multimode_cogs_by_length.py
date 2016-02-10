# TODO: this module has not been finished yet

# This file defines a method to determine the number of modes based on COG
# length distribution. It uses Kernel Density Estimation.

from import_proxy import *
from filedefs import *
from taxonomy import *
import numpy as np
from sklearn import mixture
from scipy import stats
import matplotlib.pyplot as plt

def modeCount(input, thresholdMax, thresholdValley, minPoints):
    """
    :param input: input list
    :param thresholdMax: Ignore maximums that are less than this fraction of
        the absolute max
    :param thresholdValley: valley should be less than this as portion of
        the minimum max on either side
    :param minPoints: consider multimode only if we have more points than this
    :return: number of modes
    """

    input = np.array(sorted(input))

    gkde = stats.gaussian_kde(input, bw_method='silverman')
    xCoord = np.arange(input[0], input[-1], 1.)
    yCoord = gkde.evaluate(xCoord)

    plt.plot(xCoord, yCoord)
    plt.show()
    return

    maxRatio = 0.
    bestModeCount = 1
    if (len(input) / minPointsPerCluster < maxClusters):
        maxClusters = len(input) / minPointsPerCluster
    for n in range(2, maxClusters + 1):
        ratio = clusterizeOneDimension(input, n, threshold)
        if ratio > maxRatio:
            maxRatio = ratio
            if maxRatio >= threshold:
                bestModeCount = n

    return bestModeCount


def clusterizeOneDimension(input, clusterCount, thresh):
    inputs = input.reshape([-1,1])
    g = mixture.GMM(n_components=clusterCount)
    g.fit(inputs)
    if not g.converged_:
        return -1.0
    logProb, resp = g.score_samples(inputs)
    print "Sum of prob: ", sum(logProb)
    print "Means: ", g.means_
    print "Weights: ", g.weights_
    return (float(len(filter(lambda x: x > thresh, [max(x) for x in resp]))) /
        len(inputs))


#Build histogram
if __name__ == "__main__":

    print("reading COG instance set...")
    with open(SAMPLE_COG_INST_SET(), 'r') as fset:
        cogInstSet = json.load(fset, object_hook = UtilJSONDecoderDictToObj)
    print("Read %d COG instances" % len(cogInstSet))

    cogLenDict = {}
    for cogInst in cogInstSet:
        l = cogLenDict.get(cogInst.getName(), [])
        l.append(cogInst.getLen())
        cogLenDict[cogInst.getName()] = l

    for cogName, lenList in cogLenDict.items():
        print("Processing %s" % cogName)
        modeCount(lenList, 0.95, 3, 6)








