# This module builds cleaned prokaryot dictionary: directory -> ProkDnaSet

from import_proxy import *
from taxonomy import TaxonomyParser, Taxa
from filedefs import *

cleanDict = {}

with open(PROK_GENOME_DICT(), 'r') as fdict:
    masterDict = json.load(fdict, object_hook = UtilJSONDecoderDictToObj)

for d,pg in masterDict.iteritems():
    sn = pg.getStrainList()
    if len(sn) == 1:
        cleanDict[d] = pg.getStrain(sn[0])
    else:
        print("%d strains in %s" % (len(sn), d))

with open(PROK_CLEAN_GENOME_DICT(), 'w') as fdict:
    json.dump(cleanDict, fdict, cls = UtilJSONEncoder, sort_keys = True,
              indent = 4)

print("%s: output %d entries" % (PROK_CLEAN_GENOME_DICT(), len(cleanDict)))

# Now creating files with Taxonomy
taxonomyParser = TaxonomyParser(config.TAXONOMY_FILE())

for pds in cleanDict.values():
    taxonomyParser.addProkDnaSet(pds)

taxonomyParser.process()
print(taxonomyParser.stats())