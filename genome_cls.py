# This file defines classes used in the Prokaryotic Phylogeny Project.

import re
from import_proxy import *

class ProkDna(UtilObject):
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
    # "complete chromsome" needs to be removed, otherwise it consumes
    # the next word as the name of this chromosome
    patCompleteChromosome = re.compile(r'complete chromosome ')

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

    def processPattern(self, name, pattern, default):
        # Extracts matched feature, and returns updated DNA name
        # excluding the pattern, and teh feature value (or default)
        m = pattern.match(name)
        if not m:
            return (name, default)
        featureName = m.group(3)
        if not featureName:
            featureName = "NONAME"
        # Exclude the found feature from the name, and trim white space
        name = ProkDna.patWhiteSpace.sub(' ', m.group(1) + m.group(4))
        return (ProkDna.patWhiteSpaceComma.sub(',', name), featureName)

    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self.fullPttName = kwargs["fullPttName"]
        _, self.ptt = os.path.split(self.fullPttName)

        with open(self.fullPttName, 'r') as f:
            # Take the first line
            for l in f:
                # Get DNA name
                self.name = ProkDna.removeQuotes(l.strip().lower())
                self.name = ProkDna.patWhiteSpace.sub(' ', self.name)
                self.name = ProkDna.patCompleteChromosome.sub('', self.name)
                break

        self.chr = None
        self.strain = None
        self.isClone = False
        self.phage = None
        self.plasmid = None

        if ProkDna.patClone.match(self.name):
            self.isClone = True
            return

        self.name, self.phage = self.processPattern(self.name,
                                                    ProkDna.patPhage, None)
        self.name, self.plasmid = self.processPattern(self.name,
                                                  ProkDna.patPlasmid, None)
        if self.plasmid or self.phage:
            return

        self.name, self.chr = self.processPattern(self.name,
                                               ProkDna.patChromosome, "1")
        self.xlateChromosomeStr()
        # Some names have an extra chromosome string, remove it
        self.name, _ = self.processPattern(self.name, ProkDna.patChromosome,
                                           None)

        defaultStrainName = "MAIN"
        self.name, self.strain = self.processPattern(self.name,
                                                    ProkDna.patStrain1,
                                                    defaultStrainName)
        if self.strain == defaultStrainName:
            # Try another pattern
            self.name, self.strain = self.processPattern(self.name,
                                                        ProkDna.patStrain2,
                                                        defaultStrainName)

        # Leave only the name up to the comma
        self.name = ProkDna.patWhiteSpaceComma.sub(',', self.name)
        self.name = self.name.split(',', 1)[0]


    def getPtt(self):
        return self.ptt

    def getPttBase(self):
        return self.ptt.rpartition('.')[0]

    def getFullPttName(self):
        return self.fullPttName

    def getDir(self):
        return self.fullPttName.split('/')[-2]

    def key(self):
        return self.getPttBase() + "-" + self.getDir()

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


class ProkDnaSet(UtilObject):
    """
    Collection of ProkDna objects, indexed by chromosome id (0 based)
    Attributes:
        strain - strain of the main DNA
        dct - a dictionary of chromosome id -> ProkDna mappings
        name - offcial name of this strain, from the first line of the PTT
            file
    """

    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self.strain = None
        self.name = None
        self.dct = {}
        self.dir = None

    def add(self, prokDna):
        strain = prokDna.getStrain()
        if not self.strain:
            self.strain = strain
        if self.strain != strain:
            raise UtilError("ProkDnaSet %s and ProkDna %s got strain "
                              "mismatch" % (self, prokDna))
        if not self.name:
            self.name = prokDna.getName()
        if self.name != prokDna.getName():
            raise UtilError("Name mismatch in ProkDnaSet: %s vs %s, ProkDna "
                            "file %s" %
                            (self.name, prokDna.getName(),
                             prokDna.getFullPttName()))

        if not self.dir:
            self.dir = prokDna.getDir()
        if self.dir != prokDna.getDir():
            raise UtilError("Directory mismatch in ProkDnaSet: %s vs %s, "
                            "ProkDna file %s" %
                            (self.dir, prokDna.getDir(),
                             prokDna.getFullPttName()))

        chromId = prokDna.getChromId()
        if chromId in self.dct:
            raise UtilError("ProkDnaSet %s adds ProkDna %s with the same "
                              "chromosome" % (self, prokDna))
        self.dct[chromId] = prokDna

    def getName(self):
        return self.name

    def getDir(self):
        return self.dir

    def getStrain(self):
        return self.strain

    def key(self):
        return self.strain + "-" + self.dir

    def getChromCount(self):
        return len(self.dct)

    def getChrom(self, id):
        return self.dct[id]

    def getChromIdList(self):
        return self.dct.keys()


class ProkGenome(UtilObject):
    """
    Describes full organism's genome
    Attributes:
        dir - directory of this genome
        chroms - DNA from chromosomes, a dictionary of <strain> -> ProkDnaSet
        phages - list of ProkDna corresponding to phages
        plasmids - list of ProkDna corresponding to plasmids
        clones - list of clones
    """

    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self.dir = kwargs["dir"]
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
                raise UtilError("ProkDnaSet %s: inconsistent chromosome "
                                  "count" % self.chroms)
            if count > 1:
                if not chromIdSet:
                    chromIdSet = set(prokDnaSet.getChromIdList())
                else:
                    if chromIdSet != set(prokDnaSet.getChromIdList()):
                        raise UtilError("ProkDnaSet %s: inconsistent "
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


class ProkCog(UtilObject):
    """
    Class describing a prokaryote COG.
    Attributes:
        name - name of the COG
        chrom - ProkDna.key for this COG instance
        pttLine - line in the PTT file
        strand - strand of the chromosome
        start - strating position in the chromosome
        len - length, in terms of proteins
        cntLine - chain of aminiacids in the protein, as a line number in
            the <COGNAME>.cnt file in the work files directory
    """

    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self.__dict__.update(kwargs)

    def key(self):
        return str(self.pttLine) + "-" + self.chrom

    def getName(self):
        return self.name

    def getChrom(self):
        return self.chrom





