# This module builds set of CogInst, and Cog objects

import re
from import_proxy import *
from filedefs import *

cogPat = re.compile(r'^COG.*')

# Set of CogInst
fullCogInstSet = set()

# Dictionary mapping COG name -> Cog
fullCogDict = {}

# Dictionary of FASTA files containing COG proteins, file name -> current
# line number
faFileDict = {}

# Dictionary of genome -> COG count
genomeDict = {}

# Use multifile to append files
multiFile = UtilMultiFile(3000, "a")

# Files with COG protein sequences
def cogFastaFileName(nameBase):
    return config.WORK_FILES_DIR() + nameBase + ".fa"

# Returns a set of COGs contaned in this prokDna
def getCogSet(prokDna):
    # Dictionary of PID -> protein
    cogProteinDict = {}

    # Read all COG proteins from the corresponding faa file
    protein = ""
    pid = None
    faaFileName = prokDna.getFullPttName().rpartition('.')[0] + ".faa"
    with open(faaFileName, 'r') as ffaa:
        for l in ffaa:
            l = l.strip()
            ll = l.split('|')
            if len(ll) == 5:
                # Protein descriptor
                if pid:
                    cogProteinDict[pid] = protein
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

    if pid:
        cogProteinDict[pid] = protein

    return buildCogSet(prokDna, cogProteinDict)

def buildCogSet(prokDna, cogProteinDict):
    cogInstSet = set()

    with open(prokDna.getFullPttName(), 'r') as fptt:
        # Skip first 3 lines
        lineCount = 0
        for _ in fptt:
            lineCount += 1
            if lineCount == 3:
                break

        for l in fptt:
            lineCount += 1
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
                    prokDna.getFullPttName(), lineCount))
                continue

            if cogPid not in cogProteinDict:
                print("COG from file %s line %u pid %s not in FAA file" % (
                    prokDna.getFullPttName(), lineCount, cogPid))
                continue

            faFileName = cogFastaFileName(cogName)
            faLineNumber = faFileDict.get(faFileName, 1)

            cogInst = CogInst(name = cogName, chrom = prokDna.key(), pttLine =
                lineCount, strand = cogStrand, start = cogStart, len = cogLen,
                faLine = faLineNumber)
            cogInstSet.add(cogInst)

            multiFile.write(faFileName, '>' + cogInst.key() + '\n' +
                            cogProteinDict[cogPid] + '\n')
            faFileDict[faFileName] = faLineNumber + 2

    return cogInstSet


with open(PROK_CLEAN_GENOME_DICT(), 'r') as fdict:
    masterDict = json.load(fdict, object_hook = UtilJSONDecoderDictToObj)

for d, prokDnaSet in masterDict.iteritems():
    for cid in prokDnaSet.getChromIdList():
        prokDna = prokDnaSet.getChrom(cid)
        cogInstSet = getCogSet(prokDna)
        if cogInstSet:
            fullCogInstSet.update(cogInstSet)

multiFile.closeAll()
print ("MultiFile stat: %s" % multiFile.getStats())

# Now build fullCogDict
for cogInst in fullCogInstSet:
    cog = fullCogDict.get(cogInst.getName(), Cog(name=cogInst.getName()))
    cog.addCogInst(cogInst)
    fullCogDict[cogInst.getName()] = cog

cogNamesList = fullCogDict.items()
for name, cog in cogNamesList:
    genomes = cog.calculate()
    assert(cog.getGenCount() >= 1)
    for g in genomes:
        genomeDict[g] = genomeDict.get(g, 0) + 1

print("%d Cog Instances" % len(fullCogInstSet))
print("Dumping to a file...")
with open(COG_INST_SET(), 'w') as f:
    json.dump(fullCogInstSet, f, cls = UtilJSONEncoder, sort_keys = True,
              indent = 4)

print("%d Cogs total" % len(fullCogDict))
print("Dumping to a file...")
with open(COG_LIST(), 'w') as f:
    cogList = sorted(fullCogDict.values(), key = lambda x: x.instCount,
                     reverse = True)
    json.dump(cogList, f, cls = UtilJSONEncoder, sort_keys=True, indent=4)

print("%d genomes got COGs" % len(genomeDict))
with open(GENOME_COG_CNT_LIST(), 'w') as f:
    json.dump(sorted(genomeDict.items(), key = lambda x: x[1]), f, cls =
        UtilJSONEncoder, sort_keys=True, indent=4)

# See how many genomes got both COGs and taxa
with open(PROK_TAXA_DICT(), 'r') as f:
    taxaDict = json.load(f, object_hook = UtilJSONDecoderDictToObj)
taxaSet = set(taxaDict.keys())
genomeWithCogsSet = set(genomeDict.keys())
print("%d genomes with both taxa and COGs" % len(set.intersection(taxaSet,
    genomeWithCogsSet)))


