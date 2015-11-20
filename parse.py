import os
import re
import sys
import copy
import numpy
import operator
import math

bactDict = dict()
cogsDict = dict()
cogsLengths = dict()
parse_error_count = 0
totalLinesProcessed = 0

ptt_pat = re.compile(r'.*\.ptt')

def main():
	for dirName, subdirList, fileList in os.walk("../BACTERIA"):
		if dirName in [".", "./PROC", "./FILES"]:
			continue        
		d = dict()
		bactDict[dirName] = d
		for fname in fileList:
			if ptt_pat.match(fname):
				processCogs(dirName + "/" + fname, d)
	temp = bactDict.keys()
	for d in temp:
		if len(bactDict[d]) == 0:
			del bactDict[d]
	reverseDict()
	buildCogsLengths()	
	processResults()

def buildCogsLengths():
	print("Building cogsLenghths")
	for c, bactDict in cogsDict:
		l = []
		for b, entries in bactDict:
			l += [e[0] for e in entries]
		l.sort()
		cogsLengths[c] = l

def reverseDict():
	print("Building cogs -> [bacts]")
	for b, cogs in bactDict.iteritems():
		for cog, cogEntries in cogs.iteritems():
			bactPerCog = cogsDict.get(cog, {})
			bactPerCog[b] = cogEntries 
			cogsDict[cog] = bactPerCog
			

def processResults():
	cogsPerBactList = []
	bactsPerCogList = []
	cogsInOneBactLogList = []
	for b, cogs in bactDict.iteritems():
		cogsPerBactList.append(math.log(len(cogs), 1.5))
		for c, entries in cogs.iteritems():
			cogsInOneBactLogList.append(len(entries).bit_length())
	for c, bacts in cogsDict.iteritems():
		bactsPerCogList.append(math.log(len(bacts), 1.5))

	hist = numpy.histogram(cogsInOneBactLogList, bins = range(20))
	print("Cogs In One Bact Log Hist:\n%s\n\n" % \
		'  '.join(("%d: %d" % (indx, p)) \
		for (indx, p) in enumerate(hist[0].tolist()))) 

	hist = numpy.histogram(cogsPerBactList, bins = range(30))
	print("Cogs Per Bact HistL\n")
	for indx, v in enumerate(hist[0].tolist()):
		'  '.join(("%d: %d-%d %d" % (indx, b0, b1, p)) \
		for (indx, bp)


def processCogs(fname, cogsPerBact):
	#print("Processing %s" %  fname)
	with open(fname, 'r') as f:
		lineCount = 0
		for l in f:
			lineCount += 1
			if lineCount > 3:
				processCogsLine(fname, lineCount, l, cogsPerBact)

def progressIndicator():
	global totalLinesProcessed
	totalLinesProcessed += 1
	if totalLinesProcessed % 100000 == 0:
		print("Lines processed %d" % totalLinesProcessed)

def processCogsLine(fname, lineCount, l, cogsPerBact):
	progressIndicator()
	global parse_error_count
	ll = l.split('\t')
	if len(ll) >= 6:
		try:
			cogLen = int(ll[2])
		except:
			print("Parse error file %s line %d" % (fname, lineCount))
			parse_error_count += 1
			return
		if cogLen == 0:
			return
		cogName = ll[7]
		if cogName == "-":
			return
		entries = cogsPerBact.get(cogName, [])
		entries.append((cogLen, fname, lineCount))
		entries.sort(key=operator.itemgetter(0))
		cogsPerBact[cogName] = entries

if __name__ == "__main__":
	main()
