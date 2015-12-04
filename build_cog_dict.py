# This module builds set of ProkCogs

import re
from import_proxy import *
from filedefs import *

cogPat = re.compile(r'^COG.*')

# Set of ProkCogs
fullCogSet = set()

# Dictionary of FASTA files containing COG proteins, file name -> current
# line number
faFileDict = {}

multiFile = UtilMultiFile(3000, "a")

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
    cogSet = set()

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

            faFileName = config.WORK_FILES_DIR() + cogName + ".fa"
            faLineNumber = faFileDict.get(faFileName, 1)

            cog = ProkCog(name = cogName, chrom = prokDna.key(), pttLine =
                lineCount, strand = cogStrand, start = cogStart, len = cogLen,
                faLine = faLineNumber)
            cogSet.add(cog)

            multiFile.write(faFileName, '>' + cog.key() + '\n' +
                            cogProteinDict[cogPid] + '\n')
            faFileDict[faFileName] = faLineNumber + 2

    return cogSet


with open(PROK_CLEAN_GENOME_DICT(), 'r') as fdict:
    masterDict = json.load(fdict, object_hook = UtilJSONDecoderDictToObj)

for d, prokDnaSet in masterDict.iteritems():
    for cid in prokDnaSet.getChromIdList():
        prokDna = prokDnaSet.getChrom(cid)
        cogSet = getCogSet(prokDna)
        if cogSet:
            fullCogSet.update(cogSet)

multiFile.closeAll()
print ("MultiFile stat: %s" % multiFile.getStats())

with open(PROK_COGS_SET(), 'w') as fdict:
    json.dump(fullCogSet, fdict, cls = UtilJSONEncoder, sort_keys = True,
              indent = 4)

print("%d ProkDna's mapped to cog sets" % len(fullCogSet))
print("Created %d COG files" % len(faFileDict))
