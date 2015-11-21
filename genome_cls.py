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

class Object:
    """
    Base class defining JSON serialization methods.
    """

    def to_JSON(self, fileName):
        with open(fileName, 'w') as f:
            json.dump(self, f, default=lambda o: o.__dict__,
                      sort_keys=True, indent=4)

    @classmethod
    def from_JSON(cls, fileName):
        with open(fileName, 'w') as f:
            return type(cls.__name__, (Object,), json.load(f))

    def __repr__(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True)

    def __str__(self):
        return repr(self)



class ProkDna(Object):
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
    """

    # Legitimate chromosome representation (might have no name)
    patChromosome = re.compile(r'(.*)(\bchromosome\b(?:\s+([^\s]+))?)(.*)')
    chromosomeXlator = {"i":"1", "ii":"2", "iii":"3"}

    # Phage representation (might have no name)
    patPhage = re.compile(r'(.*)(\bphage\b(?:\s+([^\s]+))?)(.*)')

    # Strain representation (must have name)
    patStrain = re.compile(r'(.*)(\bstr[\.|ain]\b\s+([^\s]+ ))(.*)')

    # Plasmid and megaplasmid patterns (might have no name)
    patPlasmid = re.compile(r'(.*)(\b(?:mega)?plasmid\b(?:\s+([^\s]+))?)(.*)')

    patWhiteSpace = re.compile(r'\s+')
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
        return ProkDna.patWhiteSpace.sub(' ', m.group(1) + m.group(4))

    def __init__(self, pttFileName):
        self.fullPttName = pttFileName
        _, self.ptt = os.path.split(pttFileName)

        with open(pttFileName, 'r') as f:
            # Take the first line
            for l in f:
                # Get DNA name
                self.name = ProkDna.removeQuotes(l.strip().lower())
                break

        self.chr = None
        self.strain = None

        self.name = self.processPattern(self.name, ProkDna.patPhage,
                                        "phage", None)
        self.name = self.processPattern(self.name, ProkDna.patPlasmid,
                                        "plasmid", None)

        if not (self.plasmid or self.phage):
            self.name = self.processPattern(self.name, ProkDna.patChromosome,
                                            "chr", "1")
            self.xlateChromosomeStr()
            self.name = self.processPattern(self.name, ProkDna.patStrain,
                                            "strain", "main")

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


class ProkDnaSet(Object):
    """
    Collection of ProkDna objects, indexed by chromosome id (0 based)
    Attributes:
        strain - strain of the main DNA
        dict - a dictionary of chromosome id -> ProkDna mappings
    """

    def __init__(self):
        self.strain = None
        self.dict = {}

    def add(self, prokDna):
        strain = prokDna.getStrain()
        if not self.strain:
            self.strain = strain
        if self.strain != strain:
            raise GenomeError("ProkDnaSet %s and ProkDna %s got strain "
                              "mismatch" % (self, prokDna))
        chromId = prokDna.getChromId()
        if chromId in self.dict:
            raise GenomeError("ProkDnaSet %s adds ProkDna %s with the same "
                              "chromosome" % (self, prokDna))
        self.dict[chromId] = prokDna

    def getStrain(self):
        return self.strain

    def getChromCount(self):
        return len(self.dict)

    def getChrom(self, id):
        return self.dict[id]

class ProkGenome(Object):
    """
    Describes full organism's genome
    Attributes:
        dir - directory of this genome
        chroms - DNA from chromosomes, a dictionary of <strain> -> ProkDnaSet
        phages - list of ProkDna corresponding to phages
        plasmids - list of ProkDna corresponding to plasmids
    """

    def __init__(self, dir):
        self.dir = dir
        self.chroms = {}
        self.phages = []
        self.plasmids = []

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

        assert(prokDna.getPlasmid())
        self.plasmids.append(prokDna)

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

