# This file builds additional files, not needed by this project,
# but requested by other folks

from import_proxy import *
from filedefs import *

with open(PROK_TAXA_DICT(), 'r') as fdict:
    masterDict = json.load(fdict, object_hook = UtilJSONDecoderDictToObj)

nameDict = {}
for dir, taxa in masterDict.items():
    nameDict[taxa.name] = dir

dirDict = {}
names = sorted(nameDict.keys())
with open(config.WORK_FILES_DIR() + "index_name.txt", 'w') as f:
    for i, name in enumerate(names):
        # Make it 1 based
        i += 1
        f.write(str(i) + '\t' + name + '\n')
        dirDict[nameDict[name]] = i

print("reading COG instance set...")
with open(COG_INST_SET(), 'r') as fset:
    cogSet = json.load(fset, object_hook = UtilJSONDecoderDictToObj)

print("Processing cogInstList...")
cogInstList = []
dirExceptions = set()
for cogInst in cogSet:
    dir = ProkDna.getDirFromKey(cogInst.getChrom())
    if dir in dirDict:
        cogInstList.append((dirDict[dir], cogInst.getName(), cogInst.getLen()))
    else:
        if dir not in dirExceptions:
            print("Unmatched dir %s" % dir)
            dirExceptions.add(dir)

cogInstList = sorted(cogInstList, key = lambda x: x[2])
cogInstList = sorted(cogInstList, key = lambda x: x[1])
cogInstList = sorted(cogInstList, key = lambda x: x[0])

print("Dumping COG lengths to a file...")
with open(config.WORK_FILES_DIR() + "index_cog_len.txt", 'w') as f:
    for t in cogInstList:
        f.write(str(t[0]) + '\t' + t[1] + '\t' + str(t[2]) + '\n')


