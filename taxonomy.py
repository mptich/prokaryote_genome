# This module adds information about taxonomy for all organisms it can find
# it for

import re
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


