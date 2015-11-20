# This file defines classes used in the Prokaryotic Phylogeny Project.

import os
import re
import json
from utils import GenomeError

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

    patChromosome = re.compile(r'(.*)(\bchromosome\b\s+([^\s]+))(.*)')
    patPhage = re.compile(r'(.*)(\bphage\b\s+([^\s]+))(.*)')
    patStrain = re.compile(r'(.*)(\bstr[\.|ain]\b\s+([^\s]+))(.*)')
    patPlasmid = re.compile(r'(.*)(\bplasmid\b\s+([^\s]+))(.*)')
    chromosomeXlator = {"1":0, "2":1, "3":2, "i":0, "ii":1, "iii":2}

    def xlateChromosomeStr(self):
        # Translates chromosome string into 0 based id
        if self.chr in ProkDna.chromosomeXlator:
            return ProkDna.chromosomeXlator[self.chr]
        raise GenomeError("ProkDna %s got untranslatable chromosome" % self)

    def processPattern(self, name, pattern, attributeName, default):
        # Extracts matched feature, and returns updated DNA name
        # excluding the pattern)
        m = pattern.match(name)
        if not m:
            setattr(self, attributeName, default)
            return name
        setattr(self, attributeName, m.group(2))
        # Exclude the found feature from the name
        return m.group(0) + m.group(3)

    def __init__(self, pttFileName):
        self.fullPttName = pttFileName
        _, self.ptt = os.path.split(pttFileName)

        with open(pttFileName, 'r') as f:
            # Take the first line
            for l in f:
                # Get DNA name, it is terminated by a comma
                self.name = l.split(',')[0].lower()
                break

        self.chr = 0
        self.strain = None

        self.name = self.processPattern(self.name, ProkDna.patPhage,
                                        "phage", None)
        self.name = self.processPattern(self.name, ProkDna.patPlasmid,
                                        "plasmid", None)

        if not (self.plasmid or self.phage):
            self.name = self.processPattern(self.name, ProkDna.patChromosome,
                                            "chr", "1")
            self.chr = self.xlateChromosomeStr()
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

    def validate(self):
        # Make sure that all chromosomes are present
        for i in range(len(self.dict)):
            if i not in self.dict:
                raise GenomeError("ProkDnaSet %s chromsome %d is missing" %
                                  (self, i))

class ProkGenome(Object):
    """
    Describes full organism's genome
    Attributes:
        dir - directory of this genome
        main - main DNA, a dictionary of <strain> -> ProkDnaSet
        phages - list of ProkDna corresponding to phages
        plasmids - list of ProkDna corresponding to plasmids
    """

    def __init__(self, dir):
        self.dir = dir
        self.main = {}
        self.phages = []
        self.plasmids = []

    def add(self, prokDna):
        if self.dir != prokDna.getDir():
            raise GenomeError("ProkGenome %s and ProkDna %s: dir mismatch" %
                              (self, prokDna))
        if prokDna.getStrain():
            strain = prokDna.getStrain()
            pds = self.main.get(strain, ProkDnaSet())
            pds.add(prokDna)
            self.main[strain] = pds
            return

        if prokDna.getPhage():
            self.phages.append(prokDna)
            return

        assert(prokDna.getPlasmid())
        self.plasmids.append(prokDna)

    def verify(self):
        for strain, prokDnaset in self.main.iteritems():
            prokDnaset.verify()

    def getDir(self):
        return self.dir

    def getStrainList(self):
        return self.main.keys()

    def getStrain(self, strain):
        return self.main[strain]

    def getPhageCount(self):
        return len(self.phages)

    def getPhage(self, ind):
        return self.phages[ind]

    def getPlasmidCount(self):
        return len(self.plasmids)

    def getPlasmid(self, ind):
        return self.plasmids[ind]

