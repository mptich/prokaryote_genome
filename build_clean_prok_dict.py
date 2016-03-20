# This module builds cleaned prokaryote dictionary:
# directory -> ProkDnaSet

from taxonomy import TaxonomyParser
from filedefs import *
from shared.pyutils.utils import *

cleanDict = {}

masterDict = UtilLoad(PROK_GENOME_DICT())

for d,pg in masterDict.iteritems():
    sn = pg.getStrainList()
    if len(sn) == 1:
        cleanDict[d] = pg.getStrain(sn[0])
    else:
        print("%d strains in %s" % (len(sn), d))

UtilStore(cleanDict, PROK_CLEAN_GENOME_DICT())

print("%s: output %d entries" % (PROK_CLEAN_GENOME_DICT(), len(cleanDict)))

# Now creating files with Taxonomy
taxonomyParser = TaxonomyParser(config.TAXONOMY_FILE())

for pds in cleanDict.values():
    taxonomyParser.addProkDnaSet(pds)

taxonomyParser.process()
print(taxonomyParser.stats())