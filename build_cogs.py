# This module builds set of CogInst, and Cog objects

import re
from filedefs import *
import random
from shared.pyutils.utils import *
from genome_cls import ProkDna, ProkDnaSet, ProkGenome, CogInst, Cog

cogPat = re.compile(r'^COG.*')

# Set of CogInst
fullCogInstList = []

# Sample list of CogInst, for debugging
sampleCogInstList = []

# Dictionary mapping COG name -> Cog
fullCogDict = {}

# Dictionary of FASTA files containing COG proteins, file name -> current
# line number
faFileDict = {}

# Dictionary of genome -> COG count
genomeDict = {}

# Use multifile to append files
multiFile = UtilMultiFile(100, "a")

# Only proteins containing these letters are taken into consideration
validAminoAcidSet = set("ARNDCQEGHILKMFPSTWYVBZ")

# IDs of COGs which have bad protein strings
idsOfBadProteins = set()

# IDs of COGs with missing protein strings
idsOfMissingProteins = set()

# Files with COG protein sequences
def cogFastaFileName(nameBase):
    return config.WORK_FILES_DIR() + nameBase + ".fa"

# Checks protein string for validity
def checkProtein(protein):
    return (set(protein) <= validAminoAcidSet)

def addProteinToDict(cogProteinDict, pid, protein, faaFileName, lineno):
    if not pid:
        return
    if pid and checkProtein(protein):
        cogProteinDict[pid] = protein
    else:
        print("File %s line %d: ignoring protein %s id %s" % (
            faaFileName, lineno, protein, pid))
        idsOfBadProteins.add(pid)

# Returns a set of COGs contaned in this prokDna
def getCogSet(prokDna):
    # Dictionary of PID -> protein
    cogProteinDict = {}

    # Read all COG proteins from the corresponding faa file
    protein = ""
    pid = None
    faaFileName = prokDna.getFullPttName().rpartition('.')[0] + ".faa"
    with open(faaFileName, 'r') as ffaa:
        for lineno, l in enumerate(ffaa, start = 1):
            l = l.strip()
            ll = l.split('|')
            if len(ll) >= 5:
                addProteinToDict(cogProteinDict, pid, protein, faaFileName,
                                 lineno)
                pid = ll[1]
                protein = ""
                continue
            if len(ll) == 1:
                # Line with a protein
                protein += l
                continue
            # Unknown line. Reset everything till the next descriptor
            print("File %s unknown line %s" % (faaFileName, l))
            pid = None

    addProteinToDict(cogProteinDict, pid, protein, faaFileName, lineno)

    return buildCogSet(prokDna, cogProteinDict)

def buildCogSet(prokDna, cogProteinDict):
    cogInstSet = set()

    with open(prokDna.getFullPttName(), 'r') as fptt:
        # Skip first 3 lines
        for lineno, l in enumerate(fptt, start = 1):
            if lineno <= 3:
                continue

            l = l.strip()
            ll = l.split('\t')
            if len(ll) != 9:
                # Corrupted line
                continue
            if not cogPat.match(ll[7]):
                # This is not a COG
                continue
            try:
                cogStart = int(ll[0].split('.', 1)[0])
                cogLen = int(ll[2])
                # For some reason, often there are 2 COG names, separated by
                # comma. Names are the same
                cogName = ll[7].split(',')[0]
                cogStrand = ll[1]
                cogPid = ll[3]
            except:
                print("Cant't parse file %s line %u" % (
                    prokDna.getFullPttName(), lineno))
                continue

            if cogPid in idsOfBadProteins:
                # Protein has been malformatted
                continue

            if cogPid not in cogProteinDict:
                print("COG from file %s line %u pid %s not in FAA file" % (
                    prokDna.getFullPttName(), lineno, cogPid))
                idsOfMissingProteins.add(cogPid)
                continue

            faFileName = cogFastaFileName(cogName)
            faLineNumber = faFileDict.get(faFileName, 1)

            cogInst = CogInst(_name = cogName, chrom = prokDna.key(), pttLine =
                lineno, strand = cogStrand, start = cogStart, _len = cogLen,
                faLine = faLineNumber)
            cogInstSet.add(cogInst)

            multiFile.write(faFileName, '>' + cogInst.key() + '\n' +
                            cogProteinDict[cogPid] + '\n')
            faFileDict[faFileName] = faLineNumber + 2

    return cogInstSet


masterDict = UtilLoad(PROK_CLEAN_GENOME_DICT())

for d, prokDnaSet in masterDict.iteritems():
    for cid in prokDnaSet.getChromIdList():
        prokDna = prokDnaSet.getChrom(cid)
        cogInstSet = getCogSet(prokDna)
        if cogInstSet:
            fullCogInstList += list(cogInstSet)

multiFile.closeAll()
print ("MultiFile stat: %s" % multiFile.getStats())

# Now build fullCogDict
print("Building fullCogDict...")
for cogInst in fullCogInstList:
    cog = fullCogDict.get(cogInst.name, Cog(_name=cogInst.name))
    cog.addCogInst(cogInst)
    fullCogDict[cogInst.name] = cog

# Make a sample subset of COGs, for debugging
print("Building sampleCogInstList...")
sampleCogNames = set()
for cogName in fullCogDict.keys():
    if random.randrange(40) == 0:
        sampleCogNames.add(cogName)
for cogInst in fullCogInstList:
    if cogInst.name in sampleCogNames:
        sampleCogInstList.append(cogInst)

cogNamesList = fullCogDict.items()
for name, cog in cogNamesList:
    genomes = cog.calculate()
    assert(cog.getGenCount() >= 1)
    for g in genomes:
        genomeDict[g] = genomeDict.get(g, 0) + 1

print("%d Cog Instances" % len(fullCogInstList))
print("Dumping to a file...")
UtilStore(fullCogInstList, COG_INST_LIST())

print("%d Sample Cog Instances" % len(sampleCogInstList))
print("Dumping to a file...")
UtilStore(sampleCogInstList, SAMPLE_COG_INST_LIST())

print("%d Cogs total" % len(fullCogDict))
print("Dumping to a file...")
cogList = sorted(fullCogDict.values(), key = lambda x: x.instCount,
    reverse=True)
UtilStore(cogList, COG_LIST())

print("%d genomes got COGs" % len(genomeDict))
UtilStore(sorted(genomeDict.items(), key = lambda x: x[1]),
    GENOME_COG_CNT_LIST())

# See how many genomes got both COGs and taxa
taxaDict = UtilLoad(PROK_TAXA_DICT())
taxaSet = set(taxaDict.keys())
genomeWithCogsSet = set(genomeDict.keys())
print("%d genomes with both taxa and COGs" % len(set.intersection(taxaSet,
    genomeWithCogsSet)))

print("Bad proteins %d, missing proteins %d" % (len(idsOfBadProteins),
                                                len(idsOfMissingProteins)))


