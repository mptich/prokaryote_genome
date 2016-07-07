# This module builds cleaned prokaryote dictionary:
# directory -> ProkDnaSet

from taxonomy import *
from filedefs import *
from shared.pyutils.utils import *
import glob
import csv

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

print("Building %s out of %s" % (TAXONOMY_FILE(), config.TAXONOMY_DIR()))
with open(TAXONOMY_FILE(), "w") as ftax:
    csvwriter = csv.writer(ftax)
    for fname in glob.glob(config.TAXONOMY_DIR() + "*"):
        with open(fname, "r") as f:
            name = None
            csvreader = csv.reader(f)
            taxonValList = [""] * TaxaType.hierarchySize()
            for ll in csvreader:
                if len(ll) != 4:
                    raise IOError("Bad line in file %s" % fname)
                name = ll[1].lower()
                taxonType = ll[2].lower()
                id = ll[3].lower()
                taxonIndex = TaxaType.hierarchy().index(taxonType) if \
                    taxonType in TaxaType.hierarchy() else -1
                if taxonIndex >= 0:
                    taxonValList[taxonIndex] = name
            # Remove square brackets from name
            outName = ""
            for c in name:
                if c in ['[', ']']:
                    c = ''
                outName += c
            outList = [id] + [taxonValList[0]] + [outName] + taxonValList[1:]
            csvwriter.writerow(outList)

print("Building manualMatchDict...")
manualMatchDict = {}
with open(config.MANUAL_TAXA_MATCH(), 'r') as f:
    csvreader = csv.reader(f)
    for ll in csvreader:
        if len(ll) != 2:
            raise IOError("Bad manual match line %s" % str(ll))
        if ll[0] in cleanDict:
            manualMatchDict[ll[0]] = ll[1]

# Now creating files with Taxonomy
taxonomyParser = TaxonomyParser(TAXONOMY_FILE(), manualMatchDict)

for pds in cleanDict.values():
    taxonomyParser.addProkDnaSet(pds)

taxonomyParser.process()
print(taxonomyParser.stats())