# This module adds information about taxonomy for all organisms it can find
# it for

import re
import sys
import operator
from filedefs import *
from shared.pyutils.utils import *
from genome_cls import *

class TaxonomyParser(UtilObject):

    def __init__(self, fileName):
        # Mapping of name terms -> names
        self.termDict = {}
        # Mapping of names to ProkDnaSet
        self.nameDict = {}
        # Set of unmatched organism names from Taxonomy
        self.unmatchedSet = set()

        # The main result: map of ProkDnaSet keys to Taxa's
        self.taxaDict = {}

        # Mapping of taxonomy names to Taxa
        self.taxaNamesDict = {}
        with open(fileName, 'r') as f:
            for l in f:
                ll = l.strip().lower().split(',')
                if len(ll) != 7:
                    continue
                ll = [x.strip() for x in ll]
                taxa = Taxa(_name = ll[2], _type = TaxaType(domain=ll[1], \
                    phylum=ll[3], cls=ll[4], order = ll[5], family=ll[6]))
                self.taxaNamesDict[ll[2]] = taxa

    def addProkDnaSet(self, prokDnaSet):
        name = prokDnaSet.name
        self.nameDict[name] = prokDnaSet
        # Break the name into the terms
        nameList = re.split(' |,', name)
        for n in nameList:
            nameSet = self.termDict.get(n, set())
            nameSet.add(UtilDict(terms = nameList, name = name))
            self.termDict[n] = nameSet

    def process(self):
        # We will match organisms by the best name match
        for orgName, taxa in self.taxaNamesDict.items():
            # Remove strain if it is present
            orgName = ProkDna.removeStrain(orgName)
            termList = [x for x in orgName.split(' ') if x != ""]
            # If the last term looks like a chromosome number - remove it
            if termList and \
                    (ProkDna.chromosomeStrToNumber(termList[-1]) != -1):
                termList = termList[:-1]
            if len(termList) < 2:
                print("Taxonomy: organism name %s is too short" % orgName)
                continue
            # Set of sets of names
            listOfSets = []
            for t in termList:
                listOfSets.append(self.termDict.get(t, set()))
            # Now we need to find intersection of all sets of prokDna name
            # lists
            commonNames = set.intersection(*listOfSets)
            if not commonNames:
                self.unmatchedSet.add(orgName)
            else:
                # Find the shortest name
                shortestNameLen = 1000 # Set it to a very high value
                shortestNameDictList = []
                for d in commonNames:
                    if len(d["terms"]) < shortestNameLen:
                        shortestNameLen = len(d["terms"])
                        shortestNameDictList = [d]
                    elif len(d["terms"]) == shortestNameLen:
                        shortestNameDictList.append(d)
                assert(shortestNameDictList)

                for d in shortestNameDictList:
                    prokDnaSet = self.nameDict[d["name"]]
                    self.taxaDict[prokDnaSet.dir] = taxa

        # Dump taxa dixionary and unmatched Taxa organism names
        UtilStore(self.taxaDict, PROK_TAXA_DICT())
        UtilStore(self.unmatchedSet, UNMATCHED_TAXA_SET())
        # Dump procDna directories that has not matched anything in taxonomy
        prokDnaDirSet = set([x.dir for x in self.nameDict.values()])
        matchedProkDnaDirSet = set(self.taxaDict.keys())
        UtilStore(prokDnaDirSet - matchedProkDnaDirSet,
            UNMATCHED_PROC_DNA_SET())

    def stats(self):
        return ("Taxonomy size %d, out of them unmatched %d; "
            "%d gemomes matched taxonomy" % (len(self.taxaNamesDict),
            len(self.unmatchedSet), len(self.taxaDict)))



class TaxaType(UtilObject):
    """
    This class describes a taxonomy classification of a
        Prokaryotic organism
    Attributes:
        domain, phylum, cls, order, family
    """

    def __init__(self, **kwargs):
        self.domain = ""
        self.phylum = ""
        self.cls = ""
        self.order = ""
        self.family = ""
        if self.buildFromDict(kwargs):
            return
        self.__dict__.update(kwargs)

    @property
    def key(self):
        return self.family + "-" + self.order + "-" + self.cls + "-" + \
            self.phylum + "-" + self.domain

    def distance(self, other):
        if self.domain != other.domain:
            return 5
        elif self.phylum != other.phylum:
            return 4
        elif self.cls != other.cls:
            return 3
        elif self.order != other.order:
            return 2
        elif self.family != other.family:
            return 1
        else:
            return 0

    @staticmethod
    def hierarchy():
        return ["domain", "phylum", "cls", "order", "family"]

    def __repr__(self):
        s = "{"
        for n in TaxaType.hierarchy():
            s += " " + n + ": " + getattr(self, n) + " "
        s += " }"
        return s

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return (self.key == other.key)

    @staticmethod
    def maxDistance():
        return 5

class Taxa(UtilObject):
    """
    This class describes a specific instance of a Prokaryotic organism
    along with its taxonomy type
    Attributes:
        _name - name of the organism
        _type - TaxaType for this organism
    """

    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self.__dict__.update(kwargs)

    @property
    def key(self):
        return self._type.key + "_" + self._name

    def distance(self, other):
        return self._type.distance(other.type)

    def __repr__(self):
        s = "< " + self._name + " " + repr(self._type) + ">"
        return s

    @property
    def type(self):
        return self._type

    @property
    def name(self):
        return self._name


class TaxaTypeNode(UtilObject):
    def __init__(self, type, name, parent):
        self.dict = {}
        self.type = type
        self.name = name
        self.parentId = id(parent)
        self.parent = parent

    @property
    def key(self):
        return (self.name, self.parentId)

    def __eq__(self, other):
        return ((self.name == other.name) and \
               (self.parentId == other.parentId))

    def __repr__(self):
        return "{ " + self.type + ": " + self.name + " parent: " + \
            (self.parent.name if self.parent else "") + " }"

    def __str__(self):
        return self.__repr__()


class TaxaTypeTree(UtilObject):
    """
    This class represents a tree of TaxaTypeNodes. Basically, it is just
    a dictionary of top level nodes.
    """
    def __init__(self, taxaTypeSet):
        self.nodeCostDict = None
        self.dict = {}
        nameList = TaxaType.hierarchy()
        for taxaType in taxaTypeSet:
            currDict = self.dict
            parentObj = None
            for name in nameList:
                nameVal = getattr(taxaType, name, None)
                if nameVal is None:
                    break
                if nameVal not in currDict:
                    currDict[nameVal] = TaxaTypeNode(
                        name, nameVal, parentObj)
                parentObj = currDict[nameVal]
                currDict = parentObj.dict

    @staticmethod
    def taxaTypeFromNodeList(nodeList):
        return TaxaType(**dict(zip(TaxaType.hierarchy(),
            [x.name for x in nodeList])))

    def bldCostDict(self, taxaTypeCostDict):
        nodeList = []
        nodeCostDict = {}
        TaxaTypeTree.recurse(self.dict, nodeList, nodeCostDict,
                taxaTypeCostDict)
        return nodeCostDict

    def optimal(self, nodeCostDict):
        """
        :param nodeCostDict:
        :return: optimal TaxaType for the given nodeCostDict
        """
        cost, bestNodeList = TaxaTypeTree.walkDown(self.dict.values(),
            nodeCostDict, [], sys.float_info.max, 100.)
        taxaType = TaxaTypeTree.taxaTypeFromNodeList(bestNodeList)
        return (taxaType, cost)

    @staticmethod
    def walkDown(nodeList, nodeCostDict, currNodeList, bestCost, alpha):

        nodeList = sorted(nodeList, key = lambda x: nodeCostDict[x].mean)
        maxCost = nodeCostDict[nodeList[0]].mean * (1.0 + alpha)
        bestNodeList = None

        for node in nodeList:
            if nodeCostDict[node].mean > maxCost:
                break
            currNodeList.append(node)
            if not node.dict:
                cost = nodeCostDict[node].mean
                if cost < bestCost:
                    bestCost = cost
                    bestNodeList = currNodeList[:]
            else:
                bestCost, candidateNodeList = TaxaTypeTree.walkDown(
                    node.dict.values(), nodeCostDict, currNodeList,
                    bestCost, alpha)
                if candidateNodeList:
                    bestNodeList = candidateNodeList
            del currNodeList[-1]

        return (bestCost, bestNodeList)


    @staticmethod
    def recurse(nodeDict, nodeList, nodeCostDict, taxaTypeCostDict):
        for node in nodeDict.values():
            nodeList.append(node)
            if not node.dict:
                taxaType = TaxaTypeTree.taxaTypeFromNodeList(nodeList)
                nodeCostDict[node] = taxaTypeCostDict[repr(taxaType)]
            else:
                TaxaTypeTree.recurse(node.dict, nodeList, nodeCostDict,
                    taxaTypeCostDict)
                nodeCostDict[node] = reduce(lambda x,y: x.combine(y),
                    [nodeCostDict[x] for x in node.dict.values()])
            del nodeList[-1]

    def utilJsonDump(self, currDict = None, nodeAttribDict = None):
        if currDict is None:
            currDict = self.dict
        localList = []
        for name, node in currDict.iteritems():
            if nodeAttribDict is None:
                localList.append(name)
            else:
                localList.append((name, nodeAttribDict[node]))
            if node.dict:
                localList.append(self.utilJsonDump(currDict = node.dict,
                    nodeAttribDict = nodeAttribDict))
        return localList

    @classmethod
    def utilJsonLoad(dumpStr):
        """
        Not implemented yet
        :return: object of this class
        """
        raise Exception("Not implemented")
        return None


