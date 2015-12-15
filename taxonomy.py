# This module adds information about taxonomy for all organisms it can find
# it for

import re
from import_proxy import *
from filedefs import *

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
                taxa = Taxa(domain=ll[1], phylum=ll[3], cls=ll[4], \
                        name = ll[2], order = ll[5], family=ll[6])
                self.taxaNamesDict[ll[2]] = taxa

    def addProkDnaSet(self, prokDnaSet):
        name = prokDnaSet.getName()
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
                    self.taxaDict[prokDnaSet.getDir()] = taxa

        # Dump taxa dixionary and unmatched Taxa organism names
        with open(PROK_TAXA_DICT(), 'w') as f:
            json.dump(self.taxaDict, f, cls = UtilJSONEncoder,
                      sort_keys = True, indent = 4)
        with open(UNMATCHED_TAXA_SET(), 'w') as f:
            json.dump(self.unmatchedSet, f, cls = UtilJSONEncoder,
                      sort_keys = True, indent = 4)

    def stats(self):
        return ("Taxonomy size %d, out of them unmatched %d; "
            "%d gemomes matched taxonomy" % (len(self.taxaNamesDict),
            len(self.unmatchedSet), len(self.taxaDict)))



class Taxa(UtilObject):
    """
    This class describes a taxa of Prokaryotic organism
    Attributes:
        domain, phylum, cls, order, family
    """

    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self.__dict__.update(kwargs)

    def key(self):
        return self.family + "-" + self.order + "-" + self.cls + "-" + \
            self.phylum + "-" + self.domain

    def distance(self, other):
        if self.family == other.family:
            return 0
        elif self.order == other.order:
            return 1
        elif self.cls == other.cls:
            return 2
        elif self.phylum == other.phylum:
            return 3
        elif self.domain == other.domain:
            return 4
        else:
            return 5


