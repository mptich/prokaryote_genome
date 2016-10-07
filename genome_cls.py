# This file defines classes used in the Prokaryotic Phylogeny Project.

import re
import numpy
import os
from shared.pyutils.utils import *
from shared.pyutils.bioutils import *

class ProkDna(UtilObject):
    """
    Prokaryote DNA Molucule.
    Attributes:
    ptt: PTT file name
    fullPttName: fully qualified PTT file name
    _name: name of the DNA described in the PTT file
    chr: chromosome index, 0 based (always 0 for phages and plasmids)
    strain: name of the strain (main if no strain name; None for phages and
        plasmids)
    phage: name of the phage (None if not phage)
    plasmid: name of the plasmid (None if not plasmid)
    isClone - bool, indicating if this is some kind of clone
    """


    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self.fullPttName = kwargs["fullPttName"]
        _, self.ptt = os.path.split(self.fullPttName)

        with open(self.fullPttName, 'r') as f:
            # Take the first line
            for l in f:
                nameParser = ProkDnaNameParser(l)
                # Take just the first line
                break

        self.chr = nameParser.chr
        self.strain = nameParser.strain
        self.isClone = nameParser.isClone
        self.isElement = nameParser.isElement
        self.phage = nameParser.phage
        self.plasmid = nameParser.plasmid
        self._name = nameParser.name

    def getPtt(self):
        return self.ptt

    def isAuxiliary(self):
        return (self.isClone or self.isElement or self.phage or self.plasmid)

    def getPttBase(self):
        return self.ptt.rpartition('.')[0]

    def getFullPttName(self):
        return self.fullPttName

    @property
    def dir(self):
        return self.fullPttName.split('/')[-2]

    @property
    def key(self):
        return self.getPttBase() + ":" + self.dir

    @staticmethod
    def getDirFromKey(key):
        keyComponents = key.split(':')
        assert(len(keyComponents) == 2)
        return keyComponents[1]

    @property
    def name(self):
        return self._name

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
        _name - offcial name of this strain, from the first line of the PTT
            file
    """

    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self.strain = None
        self._name = None
        self.dct = {}
        self._dir = None

    def add(self, prokDna):
        strain = prokDna.getStrain()
        if not self.strain:
            self.strain = strain
        if self.strain != strain:
            raise UtilError("ProkDnaSet %s and ProkDna %s got strain "
                              "mismatch" % (self, prokDna))
        if not self._name:
            self._name = prokDna.name
        if self._name != prokDna.name:
            raise UtilError("Name mismatch in ProkDnaSet: %s vs %s, ProkDna "
                            "file %s" %
                            (self._name, prokDna.name,
                             prokDna.getFullPttName()))

        if not self._dir:
            self._dir = prokDna.dir
        if self._dir != prokDna.dir:
            raise UtilError("Directory mismatch in ProkDnaSet: %s vs %s, "
                            "ProkDna file %s" %
                            (self._dir, prokDna.dir,
                             prokDna.getFullPttName()))

        chromId = prokDna.getChromId()
        if chromId in self.dct:
            raise UtilError("ProkDnaSet %s adds ProkDna %s with the same "
                              "chromosome" % (self, prokDna))
        self.dct[chromId] = prokDna

    @property
    def name(self):
        return self._name

    @property
    def dir(self):
        return self._dir

    def getStrain(self):
        return self.strain

    @property
    def key(self):
        return self.strain + ":" + self._dir

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
        _dir - directory of this genome
        chroms - DNA from chromosomes, a dictionary of <strain> -> ProkDnaSet
        phages - list of ProkDna corresponding to phages
        plasmids - list of ProkDna corresponding to plasmids
        clones - list of clones
    """

    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self._dir = kwargs["_dir"]
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
    @property
    def dir(self):
        return self._dir

    @property
    def key(self):
        return self._dir

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


class CogInst(UtilObject):
    """
    Class describing a prokaryote COG instance.
    Attributes:
        _name - name of the COG instance
        chrom - ProkDna.key for this COG instance
        pttLine - line in the PTT file
        strand - strand of the chromosome
        start - strating position in the chromosome
        _len - length, in terms of proteins
        faLine - chain of aminiacids in the protein, as a line number in
            the <COGNAME>.fa file in the work files directory
    """

    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self.__dict__.update(kwargs)

    @property
    def key(self):
        return str(self.pttLine) + ":" + self.chrom

    @staticmethod
    def getDirFromKey(key):
        return ProkDna.getDirFromKey(key.split(':', 1)[1])

    @property
    def dir(self):
        return ProkDna.getDirFromKey(self.chrom)

    @property
    def name(self):
        return self._name

    def getChrom(self):
        return self.chrom

    @property
    def len(self):
        return self._len


class Cog(UtilObject):
    """
    COG describing all COG instances
    Attributes:
        _name - name of the COG
        instCount - how many instances
        genCount - in how many genomes
        meanInstPerGen - float, mean of instances per genome
        stdInstPerGen - float STD (biased, from numpy) of instances per genome
        tempDict - temporary dictionary of genome -> # of CogInst
    """

    def __init__(self, **kwargs):
        if self.buildFromDict(kwargs):
            return
        self._name = kwargs["_name"]
        self.instCount = 0
        self.genCount = None
        self.meanInstPerGen = 0.
        self.stdInstPerGen = 0.
        self.tempDict = {}

    def addCogInst(self, cogInst):
        assert(self._name == cogInst.name)
        dir = CogInst.getDirFromKey(cogInst.key)
        self.tempDict[dir] = self.tempDict.get(dir, 0) + 1
        self.instCount += 1

    def calculate(self):
        """
        :return: List of genomes
        """
        self.genCount = len(self.tempDict)
        self.meanInstPerGen = numpy.mean(self.tempDict.values())
        self.stdInstPerGen = numpy.std(self.tempDict.values())
        dirNames = self.tempDict.keys()
        del self.tempDict
        return dirNames

    def getGenCount(self):
        return self.genCount

    @property
    def key(self):
        return self._name

    @property
    def name(self):
        return self._name






