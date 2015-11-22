# This file defines classes used in the Prokaryotic Phylogeny Project.

import os
import re
import json

# Make our exception more specific
class GenomeError(Exception):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class GenomeObject(object):
    """
    Base class defining serialization methods.
    """

    def __repr__(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True)

    def __str__(self):
        return repr(self)


class GenomeJSONEncoder(json.JSONEncoder):
    """
    Converts a python object, where the object is derived from GenomeObject,
    into an object that can be decoded using the GenomeJSONDecoder.
    """
    def default(self, obj):
        if isinstance(obj, GenomeObject):
            d = obj.__dict__
            d["__genomeobjecttype__"] = str(obj.__class__)
            return d
        else:
            return json.JSONEncoder.default(self, obj)


def GenomeJSONDecoderDictToObj(d):
    ourKey = "__genomeobjecttype__"
    if ourKey in d:
        moduleName, _, className = d.pop(ourKey).rpartition('.')
        inst = type(className.encode('ascii'), (GenomeObject,), dict((
            key.encode('ascii'), value) for key, value in d.items()))
    else:
        inst = d
    return inst


class ProkDna(GenomeObject):
    """
    Prokaryote DNA Molucule.
    Attributes:
    ptt: PTT file name
    fullPttName: fully qualified PTT file name
    name: name of the DNA described in the PTT file
    chr: chromosome index, 0 based (always 0 for phages and plasmids)
    strain: name of the strain (main if no strain name; None for phages and
        plasmids)
    phage: name of the phage (None if not phage)
    plasmid: name of the plasmid (None if not plasmid)
    isClone - bool, indicating if this is some kind of clone
    """

    # Legitimate chromosome representation (might have no name)
    patChromosome = re.compile(r'(.*)(\bchromosome\b(?:\s+([^\s,]+))?)(.*)')
    chromosomeXlator = {"i":"1", "ii":"2", "iii":"3"}

    # Phage representation (might have no name)
    patPhage = re.compile(r'(.*)(\bphage\b(?:\s+([^\s,]+))?)(.*)')

    # Strain representation (must have name)
    patStrain1 = re.compile(r'(.*)(\bstr\b\.\s+([^\s,]+))(.*)')
    patStrain2 = re.compile(r'(.*)(\bstrain\b\s+([^\s,]+))(.*)')

    # Plasmid and megaplasmid patterns (might have no name)
    patPlasmid = re.compile(r'(.*)(\b(?:mega)?plasmid\b(?:\s+([^\s,]+))?)(.*)')

    # Clone match
    patClone = re.compile(r'.*\bclone\b.*')

    patWhiteSpace = re.compile(r'\s+')
    patWhiteSpaceComma = re.compile(r'\s+,')
    patQuotedText = re.compile(r'([\'\"](.*?)[\'\"])')

    def xlateChromosomeStr(self):
        # Translates chromosome string if needed
        if self.chr in ProkDna.chromosomeXlator:
            self.chr = ProkDna.chromosomeXlator[self.chr]

    @staticmethod
    def removeQuotes(s):
        miter = ProkDna.patQuotedText.finditer(s)
        # Accumulated shift, because replacement string might be of
        # different length
        shift = 0
        for m in miter:
            startPos = m.start(1) + shift
            endPos = m.end(1) + shift
            replStr = m.group(2)
            replStr = ProkDna.patWhiteSpace.sub('_', replStr)
            shift += len(replStr) + startPos - endPos
            s = s[:startPos] + replStr + s[endPos:]
        return s

    def processPattern(self, name, pattern, attributeName, default):
        # Extracts matched feature, and returns updated DNA name
        # excluding the pattern)
        m = pattern.match(name)
        if not m:
            setattr(self, attributeName, default)
            return name
        featureName = m.group(3)
        if not featureName:
            featureName = "NONAME"
        setattr(self, attributeName, featureName)
        # Exclude the found feature from the name, and trim white space
        name = ProkDna.patWhiteSpace.sub(' ', m.group(1) + m.group(4))
        return ProkDna.patWhiteSpaceComma.sub(',', name)

    def __init__(self, pttFileName):
        self.fullPttName = pttFileName
        _, self.ptt = os.path.split(pttFileName)

        with open(pttFileName, 'r') as f:
            # Take the first line
            for l in f:
                # Get DNA name
                self.name = ProkDna.removeQuotes(l.strip().lower())
                self.name = ProkDna.patWhiteSpace.sub(' ', self.name)
                break

        self.chr = None
        self.strain = None
        self.isClone = False
        self.phage = None
        self.plasmid = None

        if ProkDna.patClone.match(self.name):
            self.isClone = True
            return

        self.name = self.processPattern(self.name, ProkDna.patPhage,
                                        "phage", None)
        self.name = self.processPattern(self.name, ProkDna.patPlasmid,
                                        "plasmid", None)
        if self.plasmid or self.phage:
            return

        self.name = self.processPattern(self.name, ProkDna.patChromosome,
                                        "chr", "1")
        self.xlateChromosomeStr()
        defaultStrainName = "MAIN"
        self.name = self.processPattern(self.name, ProkDna.patStrain1,
                                        "strain", defaultStrainName)
        if self.strain == defaultStrainName:
            # Try another pattern
            self.name = self.processPattern(self.name, ProkDna.patStrain2,
                                            "strain", defaultStrainName)


    def getPtt(self):
        return self.ptt

    def getFullPttName(self):
        return self.fullPttName

    def getName(self):
        return self.name

    def getChromId(self):
        return self.chr

    def getStrain(self):
        return self.strain

    def getPhage(self):
        return self.phage

    def getPlasmid(self):
        return self.plasmid

    def getIsClone(self):
        return self.isClone


class ProkDnaSet(GenomeObject):
    """
    Collection of ProkDna objects, indexed by chromosome id (0 based)
    Attributes:
        strain - strain of the main DNA
        dct - a dictionary of chromosome id -> ProkDna mappings
    """

    def __init__(self):
        self.strain = None
        self.dct = {}

    def add(self, prokDna):
        strain = prokDna.getStrain()
        if not self.strain:
            self.strain = strain
        if self.strain != strain:
            raise GenomeError("ProkDnaSet %s and ProkDna %s got strain "
                              "mismatch" % (self, prokDna))
        chromId = prokDna.getChromId()
        if chromId in self.dct:
            raise GenomeError("ProkDnaSet %s adds ProkDna %s with the same "
                              "chromosome" % (self, prokDna))
        self.dct[chromId] = prokDna

    def getStrain(self):
        return self.strain

    def getChromCount(self):
        return len(self.dct)

    def getChrom(self, id):
        return self.dct[id]

    def getChromIdList(self):
        return self.dct.keys()


class ProkGenome(GenomeObject):
    """
    Describes full organism's genome
    Attributes:
        dir - directory of this genome
        chroms - DNA from chromosomes, a dictionary of <strain> -> ProkDnaSet
        phages - list of ProkDna corresponding to phages
        plasmids - list of ProkDna corresponding to plasmids
        clones - list of clones
    """

    def __init__(self, dir):
        self.dir = dir
        self.chroms = {}
        self.phages = []
        self.plasmids = []
        self.clones = []

    def add(self, prokDna):
        if prokDna.getStrain():
            strain = prokDna.getStrain()
            pds = self.chroms.get(strain, ProkDnaSet())
            pds.add(prokDna)
            self.chroms[strain] = pds
            return

        if prokDna.getPhage():
            self.phages.append(prokDna)
            return

        if prokDna.getIsClone():
            self.clones.append(prokDna)
            return

        assert(prokDna.getPlasmid())
        self.plasmids.append(prokDna)

    def verify(self):
        # Check that every strain contais the same number of chromosomes, and
        # that the names of chromosemes are consistent
        count = 0
        chromIdSet = None
        for strain, prokDnaSet in self.chroms.iteritems():
            if count == 0:
                count = prokDnaSet.getChromCount()
            if count != prokDnaSet.getChromCount():
                raise GenomeError("ProkDnaSet %s: inconsistent chromosome "
                                  "count" % self.chroms)
            if count > 1:
                if not chromIdSet:
                    chromIdSet = set(prokDnaSet.getChromIdList())
                else:
                    if chromIdSet != set(prokDnaSet.getChromIdList()):
                        raise GenomeError("ProkDnaSet %s: inconsistent "
                                          "chromosome names" % self.chroms)

    def getDir(self):
        return self.dir

    def getStrainList(self):
        return self.chroms.keys()

    def getStrain(self, strain):
        return self.chroms[strain]

    def getPhageCount(self):
        return len(self.phages)

    def getPhage(self, ind):
        return self.phages[ind]

    def getPlasmidCount(self):
        return len(self.plasmids)

    def getPlasmid(self, ind):
        return self.plasmids[ind]

    def getCloneCount(self):
        return len(self.clones)

    def getClone(self, ind):
        return self.clones[ind]

