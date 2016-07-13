from taxonomy import *
import sys
from shared.pyutils.distance_matrix import *


def createCogDict(cogLengthFilter):

    cogLengthDict = DefDict(list)

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

    print("Storing cogDict...")
    UtilStore(cogDict, COG_DICT())

if __name__ == "__main__":

    if len(sys.argv) == 1:
        cogLengthFilter = 3.0
    else:
        cogLengthFilter = float(sys.argv[1])
    print ("Using cogLengthFilter %f" % cogLengthFilter)

    createCogDict(cogLengthFilter)

