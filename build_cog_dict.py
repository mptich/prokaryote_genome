# This module builds 2 dicionaries:
# COGNAME -> set of ProkDna's
# ProkDna -> set of ProkCog's


from import_proxy import *

cogDict = {}
chromDict = {}

def buildCogList(prok):

with open(config.PROK_CLEAN_GENOME_DICT(), 'r') as fdict:
    masterDict = json.load(fdict, object_hook = UtilJSONDecoderDictToObj)

for d, pds in masterDict:
    for cid in pds.getChromIdList():
        prokDna = pds.getChrom(cid)
        cogList = buildCogList(prokDna)
        for cog in cogList:
            pDnaSet = cogDict.get(cog.name, set())
            pDnaSet


