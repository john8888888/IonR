#!/usr/bin/env python
import os;
import sys;
import fileinput;
from collections import defaultdict;
from numpy import *;


def Set_Chr_Nr_ (CHR):
    """ Sort by chromosome """
    if chr: 
        New = CHR[3:]
        if New == 'X': New = 23
        elif New == 'Y': New = 24
        elif New == 'M': New = 25
        else: New = int(New)
    else:
        New = 0
    return New



if len(sys.argv) != 3:
	sys.exit("Usage: ion_generateStats.py <coverage.name> <merged.bed.file>");


coverageName = sys.argv[1];
mergedBedFile = sys.argv[2];

#num of base in target, now read in merged bed file and get total base
#instead of reading from targetRecords, which won't deal well with un-merged bed file
numBasesInTarget = 0;
numBasesInTargetByChr = {};
numBasesInTargetByChr = defaultdict(float)
numBasesInTargetAutoChr = 0

for line in fileinput.input(mergedBedFile):
	tokens = line.strip().split("\t");
	entryLen = (int(tokens[2]) - int(tokens[1]));
	numBasesInTarget += entryLen;
	cid = tokens[0].upper()
	if ((cid != "CHRX") & (cid != "CHRY")):
	#	print "tokens[0] added to ByChr targets"
		numBasesInTargetAutoChr += entryLen;

	numBasesInTargetByChr[cid] += entryLen;

coverageLineCt = 0;
coverageLineCtByChr = {};
coverageLineCtByChr = defaultdict(int);
coverageLineCtAutoChr = 0;
#numNoBaseBias = 0;
for line in fileinput.input(coverageName):
	coverageLineCt += 1;
	line = line.strip();
	tokens = line.split("\t");
	cid = tokens[0].upper()
	if ((cid != "CHRX") & (cid != "CHRY")):
		coverageLineCtAutoChr += 1;
	coverageLineCtByChr[cid] += 1;


#percentile coverage, using method developed by ANH with what ILMN uses
numNotCoveredBase = int(numBasesInTarget - coverageLineCt);
#coverageNeeded = [];
numNotCoveredBaseAutoChr = int(numBasesInTargetAutoChr - coverageLineCtAutoChr);
coverageNeededAutoChr = [];

numNotCoveredBaseByChr = {}

for k, v in coverageLineCtByChr.iteritems():
	numNotCoveredBaseByChr[k] = int(numBasesInTargetByChr[k] - v);


coverageNeededByChr = {};

#for i in range(numNotCoveredBase):
    #coverageNeeded.append(0);
for i in range(numNotCoveredBaseAutoChr):
    coverageNeededAutoChr.append(0);
for k in numNotCoveredBaseByChr.keys():
	for i in range(numNotCoveredBaseByChr[k]):
		if k in coverageNeededByChr.keys():
			coverageNeededByChr[k].append(0);
		else:
			coverageNeededByChr[k] = []
			coverageNeededByChr[k].append(0);



for line in fileinput.input(coverageName):
    line = line.strip();
    tokens = line.split();
    d = float(tokens[2]);
    #coverageNeeded.append(d);

    cid = tokens[0].upper()
    if ((cid != "CHRX") & (cid != "CHRY")):
	    coverageNeededAutoChr.append(d);

    if cid in coverageNeededByChr.keys():
	    coverageNeededByChr[cid].append(d);
    else:
	    coverageNeededByChr[cid] = []
	    coverageNeededByChr[cid].append(d);


coverageNeededAutoChr = array(coverageNeededAutoChr);
meanCoveragePerBaseAutoChr = mean(coverageNeededAutoChr);
numBaseGT02MeanAutoChr = sum(coverageNeededAutoChr>0.2*meanCoveragePerBaseAutoChr);
percBaseGT02MeanAutoChr = 100*(float(numBaseGT02MeanAutoChr)/float(numBasesInTargetAutoChr));
#numBaseGT01Mean = sum(coverageNeeded>0.1*meanCoveragePerBase);
#percBaseGT01Mean = 100*(float(numBaseGT01Mean)/float(numBasesInTarget));
#numBaseGT001Mean = sum(coverageNeeded>0.01*meanCoveragePerBase);
#percBaseGT001Mean = 100*(float(numBaseGT001Mean)/float(numBasesInTarget));


meanCoveragePerBaseByChr = {}
meanCoveragePerBaseByChr = defaultdict(float)
numBaseGT02MeanByChr = {}
numBaseGT02MeanByChr = defaultdict(float)
percBaseGT02MeanByChr = {}
percBaseGT02MeanByChr = defaultdict(float)

print "chromosome\tmeanCoveragePerBaseByChr\tnumBaseGT02MeanByChr\tpercBaseGT02MeanByChr"
chrList = coverageNeededByChr.keys()
#chrList.append('CHRM')
#chrList.append('CHRX')
#chrList.append('CHRY')
#print chrList
chrList.sort(lambda x,y: x-y, key=lambda x: Set_Chr_Nr_(x))
#print chrList


for k in chrList:
#for k in sorted(coverageNeededByChr.keys()):
	coverageNeededByChr[k] = array(coverageNeededByChr[k])
	meanCoveragePerBaseByChr[k] = mean(coverageNeededByChr[k])
	numBaseGT02MeanByChr[k] = sum(coverageNeededByChr[k]>0.2*meanCoveragePerBaseByChr[k])
	percBaseGT02MeanByChr[k] = 100*(float(numBaseGT02MeanByChr[k]/float(numBasesInTargetByChr[k])))
	
	print "%s\t%2.2f\t%d\t%2.2f" % (k, meanCoveragePerBaseByChr[k], numBaseGT02MeanByChr[k], percBaseGT02MeanByChr[k])
#auto chr
print "AutoChr\t%2.2f\t%d\t%2.2f" % (meanCoveragePerBaseAutoChr, numBaseGT02MeanAutoChr, percBaseGT02MeanAutoChr)


#chromosome sorting
#a = ['chr1', 'chr10', 'chr5', 'chrx']
#a.sort(lambda x,y: x-y, key=lambda x: Set_Chr_Nr_(x))
#print a

