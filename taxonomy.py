# This module adds information about taxonomy for all organisms it can find
# it for

import re
import sys
import operator
import copy
from filedefs import *
from shared.pyutils.utils import *
from shared.pyutils.UtilNormDistrib import *
from genome_cls import *
import csv

class TaxonomyParser(UtilObject):

    def __init__(self, taxaFileName, manualMatchDict):
        # Mapping of name terms -> names
        self.termDict = {}
        # Mapping of names to ProkDnaSet
        self.nameDict = {}
        # Set of unmatched organism names from Taxonomy
        self.unmatchedSet = set()

        # The main result: map of ProkDnaSet keys to Taxa's
        self.taxaDict = {}

        # Official name to dir dictionary
        self.officialNameDirDict = {}

        # Manual match count
        self.manualMatchCount = 0

        # Mapping of taxonomy names to Taxa
        self.taxaNamesDict = {}
        with open(taxaFileName, 'r') as f:
            csvreader = csv.reader(f)
            for ll in csvreader:
                if len(ll) != 11:
                    raise IOError("Bad taxonomy line %s" % str(ll))
                ll = [x.strip() for x in ll]
                taxa = Taxa(_name = ll[2], _type = TaxaType.newTaxaType(
                    ll[1], ll[3], ll[4], ll[5], ll[6], ll[7],
                    ll[8], ll[9], ll[10]))
                self.taxaNamesDict[ll[2]] = taxa

        for dir, taxaName in manualMatchDict.iteritems():
            self.taxaDict[dir] = self.taxaNamesDict[taxaName]
            self.officialNameDirDict[taxaName] = dir
            del self.taxaNamesDict[taxaName]
            self.manualMatchCount += 1

    def addProkDnaSet(self, prokDnaSet):
        if prokDnaSet.dir in self.taxaDict:
            # Has been added manually
            return
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
        for officialName, taxa in self.taxaNamesDict.items():
            # Remove strain if it is present
            nameParser = ProkDnaNameParser(officialName)
            orgName = nameParser.name
            termList = [x for x in orgName.split(' ') if x != ""]
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
                    self.officialNameDirDict[officialName] = prokDnaSet.dir

        # Dump taxa dixionary and unmatched Taxa organism names
        UtilStore(self.taxaDict, PROK_TAXA_DICT())
        UtilStore(self.officialNameDirDict, NAME_DIR_DICT())
        UtilStore(self.unmatchedSet, UNMATCHED_TAXA_SET())
        # Dump procDna directories that has not matched anything in taxonomy
        prokDnaDirSet = set([x.dir for x in self.nameDict.values()])
        matchedProkDnaDirSet = set(self.taxaDict.keys())
        UtilStore(prokDnaDirSet - matchedProkDnaDirSet,
            UNMATCHED_PROC_DNA_SET())

    def stats(self):
        return ("Taxonomy size %d, out of them unmatched %d; "
            "%d gemomes matched taxonomy, out of them %d manually" %
            (len(self.taxaNamesDict) + self.manualMatchCount,
             len(self.unmatchedSet), len(self.taxaDict),
             self.manualMatchCount))



class TaxaType(UtilObject):
    """
    This class describes a taxonomy classification of a
        Prokaryotic organism
    Attributes:
        See hierarchy()
    """

    taxonNames_ = ["superkingdom", "phylum", "class", "order", "family",
                   "genus", "species group", "species", "subspecies"]

    def __init__(self, **kwargs):
        if not self.buildFromDict(kwargs):
            self.__dict__.update(kwargs)

    def taxonValList(self):
        return [getattr(self, n) for n in TaxaType.hierarchy()]

    def taxonValListReversed(self):
        return [getattr(self, n) for n in reversed(TaxaType.hierarchy())]

    @property
    def key(self):
        return "_".join(self.taxonValList())

    def distance(self, other):
        count = TaxaType.hierarchySize()
        for x, y in zip(self.taxonValList(), other.taxonValList()):
            if x != y:
                return count
            count -= 1
        return 0

    def depth(self):
        it = (idx for (idx,val) in enumerate(self.taxonValListReversed()) \
            if val != "")
        size = TaxaType.hierarchySize()
        return size - next(it, size)

    def parent(self):
        depth = self.depth()
        if depth != 0:
            parent = copy.deepcopy(self)
            setattr(parent, TaxaType.hierarchy()[depth-1], "")
            return parent
        else:
            return self

    def isAncestor(self, anc):
        depth = anc.depth()
        if (depth > self.depth()):
            return False
        for i in range(depth):
            name = TaxaType.hierarchy()[i]
            if getattr(self, name) != getattr(anc, name):
                return False
        return True

    def commonAncestor(self, other):
        type1 = self
        type2 = other
        while type1 != type2:
            depth1 = type1.depth()
            depth2 = type2.depth()
            if depth1 >= depth2:
                type1 = type1.parent()
            if depth2 >= depth1:
                type2 = type2.parent()
        return type1

    @staticmethod
    def hierarchy():
        return TaxaType.taxonNames_

    @staticmethod
    def hierarchySize():
        return len(TaxaType.taxonNames_)

    @staticmethod
    def newTaxaType(*taxons):
        return TaxaType(**dict(zip(TaxaType.hierarchy(), taxons)))

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
        return TaxaType.hierarchySize()

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
    def __init__(self, type):
        self.children = set()
        self.dirs = set()
        self.type = type

    @property
    def key(self):
        return self.type.key

    def __eq__(self, other):
        return (self.type == other.type)

    def __repr__(self):
        return repr(self.type)

    def __str__(self):
        return self.__repr__()

    def addChild(self, type):
        self.children.add(type)

    def addDir(self, dir):
        self.dirs.add(dir)


class TaxaTypeTree(UtilObject):
    """
    This class represents a tree of TaxaTypeNodes. Basically, it is just
    a list, indexed by depth, of dictionaries mapping TaxaTypes into
    TaxaTypeNodes
    """
    def __init__(self, taxaDict):
        self.levels = []
        for i in range(TaxaType.hierarchySize() + 1):
            self.levels.append({})

        for dir, taxa in taxaDict.iteritems():
            type = taxa.type
            childNode = None
            while True:
                d = self.levels[type.depth()]
                node = d.get(type, TaxaTypeNode(type))
                d[type] = node
                node.addDir(dir)
                if childNode:
                    node.addChild(childNode)
                childNode = node
                type = type.parent()
                if type.depth() == 0:
                    break

    def getDirSet(self, type):
        return self.levels[type.depth()][type].dirs

    def getDirCount(self, type):
        return len(self.levels[type.depth()][type].dirs)

    def getChildrenSet(self, type):
        return self.levels[type.depth()][type].children

    def getNode(self, type):
        return self.levels[type.depth()][type]

    def getAllTypesSet(self):
        s = set()
        for d in self.levels:
            s |= set(d.keys())
        return s




