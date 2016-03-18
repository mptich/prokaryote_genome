# This module builds PROK_GENOME_DICT() file (see utils.py)

import glob
from filedefs import *
from shared.pyutils.utils import *
from genome_cls import ProkDna, ProkDnaSet, ProkGenome, CogInst, Cog

# Dictionary dir -> ProkGenome. Saved into PROK_GENOME_DICT file.
prokGenomeDict = {}

# Dictionary ProkDna key -> ProkDna. Saved into PROK_DNA_DICT file.
prokDnaDict = {}

inputCount = 0

def newProkGenome(dir):
    prokGenome = ProkGenome(_dir = dir)
    fullDir = config.PROKARYOTS_DIR() + dir + '/'
    pttFiles = glob.glob(fullDir + "*.ptt")
    for pttFileName in pttFiles:
        prokDna = ProkDna(fullPttName = pttFileName)
        prokDnaDict[prokDna.key()] = prokDna
        prokGenome.add(prokDna)
    prokGenome.verify()
    return prokGenome

with open(PROKARYOT_DIRS_FILE(), 'r') as fdirs:
    for dir in fdirs:
        inputCount += 1
        dir = dir.strip()
        try:
            prokGenome = newProkGenome(dir)
        except UtilError as e:
            print("UtilError: %s" % e)
            continue
        prokGenomeDict[dir] = prokGenome

with open(PROK_GENOME_DICT(), 'w') as fdict:
    json.dump(prokGenomeDict, fdict, cls = UtilJSONEncoder, sort_keys = True,
              indent = 4)

print("Input %d entries, output %d entries" % (inputCount,
                                               len(prokGenomeDict)))

with open(PROK_DNA_DICT(), 'w') as fdict:
    json.dump(prokDnaDict, fdict, cls = UtilJSONEncoder, sort_keys = True,
              indent = 4)

print("ProkDna dictionary: %d entries" % len(prokDnaDict))
