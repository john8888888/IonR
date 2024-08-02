#!/usr/bin/env python
import os;
import sys;
from matplotlib import use
use("Agg")
import matplotlib.pyplot as plt;
from scipy.interpolate import interp1d;
import locale;
import fileinput;
from collections import defaultdict;
from numpy import *;
from ampliconUtils import *;
from time import gmtime, strftime;

def parse_hsm_stats(fname):
	if not os.path.isfile(fname):
		return(0,0,0,0);

	f = open(fname, 'r');
	lines = f.readlines();
	f.close();
	percDesignCovered = lines[0].strip().split("=")[1].split()[0];
	perclocigt500reads = lines[1].strip().split("=")[1].split()[0];
	minlocicoverage = lines[2].strip().split("=")[1].split()[0];
	maxlocicoverage = lines[3].strip().split("=")[1].split()[0];
	return (percDesignCovered, perclocigt500reads, minlocicoverage, maxlocicoverage)

def ion_plot_prediction_coverage(coverage):
	numBases = float(len(coverage));
	coverage = array(sort(coverage));
	normalizedCoverage = coverage/mean(coverage);
	perc = [0.99,  0.98, 0.95];
	perc100 = [int(i * 100) for i in perc];
	xcov = [1., 20., 350.,];
	normCovVals = linspace(0, max(normalizedCoverage), 300);
	fracBasesCov = [];
	for val in normCovVals:
		fracBasesCov.append(sum(normalizedCoverage>=val)/numBases);

	ind = argsort(fracBasesCov);

	fracBases = [];
	s = "";
	vals = [];
	for i in arange(len(perc)):
		interpNormCov = interp(perc[i], array(fracBasesCov)[ind], array(normCovVals)[ind]);
		val = xcov[i]/interpNormCov;
		vals.append(val);
		s = s + "%2.0f%% at %dx requires %2.2f mean coverage\n" % (perc100[i], xcov[i], val);
    
    
	ax = plt.subplot(111);
	#plt.plot(normCovVals, fracBasesCov, color='r');
	#plt.xlabel('Normalized Coverage');
	#plt.ylabel('Fraction Bases Covered');
	#plt.title(s);
	#plt.grid();

	#ax = plt.subplot(122);
	coverage = array(coverage);
	numBases = len(coverage);
	cumsumCov = [];
	for i in range(int(max(coverage))):
		c = float(sum(coverage>=i))/float(numBases);
		cumsumCov.append(c*100);

	plt.xlabel('Depth of coverage');
	plt.ylabel('Percent of target bases');
	plt.grid();
	plt.plot(cumsumCov);
	plt.title('Overall Target Efficiency (Ion Torrent)');
	return (perc100, xcov, vals);

def conv_summary_row_to_json(instr):
	instr = instr.rstrip();
	charq = '"';
	stra = ''.join([charq, instr, charq]);
	strb = stra.lower().replace(' ', '_').replace("\t","\":\"");
	strc = ''.join(["\t", strb]);
	return ''.join([strc, ",\n"]);

def write_to_both_files(str):
	fhT.write (str);
	fhJ.write (conv_summary_row_to_json(str));

#if len(sys.argv) != 15:
#	sys.exit("Usage: ion_generateStats.py <target.stats.name> <target.file> <coverage.name> <amplicon.summary.csv> <merged.bed.file> <readLen.summary.txt> <bfmask.stats> <position_error_summary.txt> <output_dir> <output.json.results.file> <output.summary.table.file> <hsm.stats.norm.pri1> <hsm.stats.norm.pri9> <chromosome_stats.tab>");
if len(sys.argv) != 16:
	sys.exit("Usage: ion_generateStats.py <target.stats.name> <target.file> <coverage.name> <amplicon.summary.csv> <merged.bed.file> <readLen.summary.txt> <position_error_summary.txt> <output_dir> <output.json.results.file> <output.summary.table.file> <hsm.stats.norm.pri1> <hsm.stats.norm.pri9> <chromosome_stats.tab> <biased.base.txt> <read.info.file>");

targetstatsname = sys.argv[1];
targetname = sys.argv[2];
coverageName = sys.argv[3];
#representationName = sys.argv[4];
#wellFile = sys.argv[5];
#endReadFile = sys.argv[6];
ampSumFile = sys.argv[4];
#hsmPri1File = sys.argv[7];
#hsmPri9File = sys.argv[8];
mergedBedFile = sys.argv[5];
readLenSummary = sys.argv[6];
#bfmaskStats = sys.argv[7];
posErrSumFile = sys.argv[7];
ontargetdir = sys.argv[8];
outputJson = sys.argv[9];
outputTable = sys.argv[10];
hsmPri1File = sys.argv[11];
hsmPri9File = sys.argv[12];
#outputCoverage = sys.argv[14];
chrStatsFile = sys.argv[13];
biasedBaseFile = sys.argv[14];
readInfoFile = sys.argv[15];

fhJ = open(outputJson, 'w');
fhT = open(outputTable, 'w');
#fhC = open(outputCoverage, 'w');

#well.summary.file
wellNum = -1.;
wellPerc = -1.;
totalWell =  -1;
excludedWell = 0;
#if (os.path.isfile(bfmaskStats)):
#	f = open(bfmaskStats, 'r');
#	lines = f.readlines();
#	f.close();
#	for line in lines:
#		if 'Total Wells = ' in line:
#			totalWell = int(line.strip().split(" = ")[1]);
#		if 'Excluded Wells = ' in line:
#			excludedWell = int(line.strip().split(" = ")[1]);
#	wellNum = float(lines[1].strip().split(",")[1]);
#	wellPerc = float(lines[6].strip().split(",")[1]);
addWell = totalWell - excludedWell;

#HSM loci Pri1/9
percHSMPri1 = -1.;
percHSMPri9 = -1.;
if (os.path.isfile(hsmPri1File)):
	f = open(hsmPri1File, 'r');
	lines = f.readlines();
	f.close();
	#percHSMPri1 = float(lines[1].strip().split('%')[0].split('= ')[1]);
	percHSMPri1 = float(lines[0].strip());
if (os.path.isfile(hsmPri9File)):
	f = open(hsmPri9File, 'r');
	lines = f.readlines();
	f.close();
	#percHSMPri9 = float(lines[1].strip().split('%')[0].split('= ')[1]);
	percHSMPri9 = float(lines[0].strip());


#read in Chrx/y/auto stats by chr, input file generated by ion_getChromosomalBaseStats.py
#in this script, corresponding code is not removed (yet)
ChrX_meanCoveragePerBaseByChr = -1.;
ChrX_numBaseGT02MeanByChr = -1;
ChrX_percBaseGT02MeanByChr = -1.;
ChrY_meanCoveragePerBaseByChr = -1.;
ChrY_numBaseGT02MeanByChr = -1;
ChrY_percBaseGT02MeanByChr = -1.;
ChrAuto_meanCoveragePerBaseByChr = -1.;
ChrAuto_numBaseGT02MeanByChr = -1;
ChrAuto_percBaseGT02MeanByChr = -1.;
lines = []
if (os.path.isfile(chrStatsFile)):
	f = open(chrStatsFile, 'r');
	lines = f.readlines();
	f.close();
for k in lines:
	tokens = k.strip().split("\t")
	if "CHRX" in tokens[0].upper():
		ChrX_meanCoveragePerBaseByChr = float(tokens[1])
		ChrX_numBaseGT02MeanByChr = int(tokens[2])
		ChrX_percBaseGT02MeanByChr = float(tokens[3])
		#print "chrx found\t" + tokens[0]
	if "CHRY" in tokens[0].upper():
		ChrY_meanCoveragePerBaseByChr = float(tokens[1])
		ChrY_numBaseGT02MeanByChr = int(tokens[2])
		ChrY_percBaseGT02MeanByChr = float(tokens[3])
		#print "chry found\t" + tokens[0]

	if "AUTOCHR" in tokens[0].upper():
		ChrAuto_meanCoveragePerBaseByChr = float(tokens[1])
		ChrAuto_numBaseGT02MeanByChr = int(tokens[2])
		ChrAuto_percBaseGT02MeanByChr = float(tokens[3])
		#print "chrauto found\t" + tokens[0]

		
#print "time mark 1";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

#end.end.read
endReadNumfw = -1;
endReadNumrc = -1;
endReadNum = -1;
passNumOfTarget = 0; #originally >0.2mean is a pass, change it to 0.01 and added to the last column in ampliconSummary file
passNumOfTarget001 = 0; #originally >0.2mean is a pass, change it to 0.01 and added to the last column in ampliconSummary file
numIn10 = -1;
numStartOkEndEarly = -1;
#f = open(ampSumFile, 'r');
#f = open(endReadFile, 'r');
#lines = f.readlines();
#f.close();
#endReadNum = float(lines[0].strip().split("\t")[1].split('%')[0]);
#endReadNumfw = float(lines[1].strip().split("\t")[1].split('%')[0]);
#endReadNumrc = float(lines[2].strip().split("\t")[1].split('%')[0]);
for line in fileinput.input(ampSumFile):
	if not fileinput.isfirstline():
		endReadNumfw += float(line.strip().split(",")[12]);
		endReadNumrc += float(line.strip().split(",")[13]);
		endReadNum += float(line.strip().split(",")[14]);
		passNumOfTarget += int(line.strip().split(",")[11]);
		passNumOfTarget001 += int(line.strip().split(",")[27]);
		numIn10 += int(line.strip().split(",")[28]);
		numStartOkEndEarly += int(line.strip().split(",")[29]);
if (numIn10 != -1):
	numIn10 += 1; #beginning at -1; should apply this to endReadNum?
if (numStartOkEndEarly != -1):
	numStartOkEndEarly += 1;
#print "numStartOkEndEarly"
#print numStartOkEndEarly

#print "time mark 2";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

ampliconWithLessThan100Reads = [];
locale.setlocale(locale.LC_ALL, '');
targetRecords = read_bed_records(targetname);
#targetRecords use lineCt as keys now, convert to prvious format
tr = {}
for r in targetRecords:
	tr[targetRecords[r].title] = targetRecords[r]
targetRecords = tr

numTargets = len(targetRecords);
conversionRate = (float(passNumOfTarget)/float(numTargets))*100;
conversionRate001 = (float(passNumOfTarget001)/float(numTargets))*100;

#for line in fileinput.input(representationName):
for line in fileinput.input(ampSumFile):
	if not fileinput.isfirstline():
		tokens = line.strip().split(",");
		targetRecords[tokens[0]].setNumFwd(int(tokens[1]));
		targetRecords[tokens[0]].setNumRev(int(tokens[2]));
		if ((int(tokens[1]) + int(tokens[2])) <= 100):
			ampliconWithLessThan100Reads.append(tokens[0]);

#print "time mark 3";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

f = open(targetstatsname, 'r');
lines = f.readlines();
f.close();
totalReads = float(lines[0]);
mappedReads = float(lines[1]);
targetReads = float(lines[2]);
percAllOnTarget = (targetReads/totalReads)*100;

percEndReadfw = 100.*endReadNumfw/totalReads;
percEndReadrc = 100.*endReadNumrc/totalReads;
percEndRead = 100.*endReadNum/totalReads;
percReadIn10 = 100.*numIn10/totalReads;
percStartOkEndEarly = 100.*numStartOkEndEarly/totalReads;

if addWell == -1:
	percAddWell = -1.;
else:
	percAddWell = 100.*totalReads/addWell;

#print "time mark 4";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

#num of base in target, now read in merged bed file and get total base
#instead of reading from targetRecords, which won't deal well with un-merged bed file
numBasesInTarget = 0;
numBasesInTargetByChr = {};
numBasesInTargetByChr = defaultdict(float)

for line in fileinput.input(mergedBedFile):
	tokens = line.strip().split("\t");
	entryLen = (int(tokens[2]) - int(tokens[1]));
	numBasesInTarget += entryLen;
	#if tokens[0].lower() == "ByChr":
	#	print "tokens[0] added to ByChr targets"
	#	numBasesInTargetByChr += entryLen;

	numBasesInTargetByChr[tokens[0].lower()] += entryLen;

numBasesInTarget = float(numBasesInTarget);

numBasesInUnmergedTarget = 0;
for line in fileinput.input(targetname):
	tokens = line.strip().split("\t");
	if (line.find('track') == -1):
		entryLen = (int(tokens[2]) - int(tokens[1]));
		numBasesInUnmergedTarget += entryLen;
numBasesInUnmergedTarget = float(numBasesInUnmergedTarget);

#print "time mark 5";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

numNoStrandBias = 0; #3070
numNoStrandBias4060 = 0; #4060
numNoStrandBiasLess100 = 0; #3070 and <100 reads
numReads = [];
for t in targetRecords:
	rec = targetRecords[t];
	#numBasesInTarget += rec.getLength();
        den = float(rec.numfwd + rec.numrev);
        if den == 0:
            ratio = 0;
        else:	
            ratio = float(rec.numfwd) / den;

	numReads.append(rec.numfwd + rec.numrev);
	#if (ratio >= 0.3) & (ratio<=0.7):
	if (((ratio >= 0.3) and (ratio<=0.7)) or den < 10):
		numNoStrandBias +=1;
		if (den <= 100.):
			numNoStrandBiasLess100 +=1;
	#if (ratio >= 0.4) & (ratio<=0.6):
	if (((ratio >= 0.4) and (ratio<=0.6)) or den < 10):
		numNoStrandBias4060 +=1;
numReads = array(numReads);
meanNumReads = mean(numReads);
numAmpGT02Mean = sum(numReads>0.2*meanNumReads);
percAmpGT02Mean = 100*(float(numAmpGT02Mean)/float(numTargets));
numAmpGT01Mean = sum(numReads>0.1*meanNumReads);
percAmpGT01Mean = 100*(float(numAmpGT01Mean)/float(numTargets));
numAmpGT001Mean = sum(numReads>0.01*meanNumReads);
percAmpGT001Mean = 100*(float(numAmpGT001Mean)/float(numTargets));

#print "time mark 6";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

depth = [];
depthF = [];
depthR = [];
totalbasereads = 0;
numLines = 0;
maxDepth = -1;
coverageLineCt = 0;
coverageLineCtByChr = {};
coverageLineCtByChr = defaultdict(int);
numNoBaseBias = 0;
for line in fileinput.input(coverageName):
	coverageLineCt += 1;
	line = line.strip();
	tokens = line.split("\t");

	#if (tokens[0].lower() == "ByChr"):
	#	coverageLineCtByChr += 1;
	coverageLineCtByChr[tokens[0].lower()] += 1;

	d = float(tokens[2]);
	dF = float(tokens[3]);
	dR = float(tokens[4]);
	if ((dF + dR) != 0):
		dFratio = dF / (dF + dR);
		#if (dFratio >= 0.3) & (dFratio<=0.7):
		if (((dFratio >= 0.3) and (dFratio<=0.7)) or d < 10):
			numNoBaseBias +=1;
	if d > maxDepth:
		maxDepth = d;
	totalbasereads += d;
	depth.append(d);
	depthF.append(dF);
	depthR.append(dR);
	numLines +=1;

#print "time mark 7";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

percNoBaseBias = 100*(float(numNoBaseBias) / float(numLines));
depthCopy = depth;
depth = array(depth);
basesCovered = float(sum(depth>=1));
totalReads = float(lines[0]);
mappedReads = float(lines[1]);
targetReads = float(lines[2]);
percMappedOnTarget = (targetReads/mappedReads)*100;
percAllOnTarget = (targetReads/totalReads)*100;
percOffTarget = ((mappedReads - targetReads)/totalReads)*100;
percUnmapped = ((totalReads - mappedReads)/totalReads)*100;
avgBaseCovDepth = totalbasereads/numBasesInTarget;
avgBaseReadDepth = mean(depth);
sdBaseReadDepth = std(depth);
atLeast1x = (float(sum(depth>=1))/numBasesInTarget)*100;
atLeast10x = (float(sum(depth>=10))/numBasesInTarget)*100;
atLeast20x = (float(sum(depth>=20))/numBasesInTarget)*100;
atLeast100x = (float(sum(depth>=100))/numBasesInTarget)*100;
atLeast200x = (float(sum(depth>=200))/numBasesInTarget)*100;
atLeast500x = (float(sum(depth>=500))/numBasesInTarget)*100;

#print "time mark 8";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

#1x 10x 20x covarege, normalized as 100x
#mu = mean(depth);
mu = sum(depth)/numBasesInTarget; #including bases that are covered by no reads
muUnmerged = sum(depth)/numBasesInUnmergedTarget;
#print "mu and muUnmerged are %6.2f, %6.2f" % (mu, muUnmerged);
coeffNorm100 = mu/100;
depth = [];
for line in fileinput.input(coverageName):
	line = line.strip();
	tokens = line.split("\t");
	d = float(tokens[2])/coeffNorm100;
	depth.append(d);
depth = array(depth);
atLeast1xNorm100 = (float(sum(depth>=1))/numBasesInTarget)*100;
atLeast10xNorm100 = (float(sum(depth>=10))/numBasesInTarget)*100;
atLeast20xNorm100 = (float(sum(depth>=20))/numBasesInTarget)*100;
atLeast100xNorm100 = (float(sum(depth>=100))/numBasesInTarget)*100;
atLeast200xNorm100 = (float(sum(depth>=200))/numBasesInTarget)*100;
atLeast500xNorm100 = (float(sum(depth>=500))/numBasesInTarget)*100;

#percentile coverage, using method developed by ANH with what ILMN uses
numNotCoveredBase = int(numBasesInTarget - coverageLineCt);
coverageNeeded = [];
coverageNeeded = [0] * numNotCoveredBase;

#numNotCoveredBaseByChr = int(numBasesInTargetByChr - coverageLineCtByChr);
numNotCoveredBaseByChr = {}

for k, v in coverageLineCtByChr.iteritems():
	numNotCoveredBaseByChr[k] = int(numBasesInTargetByChr[k] - v);
#for k, v in numNotCoveredBaseByChr.iteritems():
#	print k, v

#print "time mark 9";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

coverageNeededByChr = {};

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
    coverageNeeded.append(d);
    cid = tokens[0].lower()
    #if (tokens[0].lower() == "ByChr"):
    #	     coverageNeededByChr.append(d);
    if cid in coverageNeededByChr.keys():
	    coverageNeededByChr[cid].append(d);
    else:
	    coverageNeededByChr[cid] = []
	    coverageNeededByChr[cid].append(d);

#print "time mark 10";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

#fig = plt.figure(figsize=(12,9));
#fig = plt.figure(figsize=(20,10));
#(percCov, xCov, estCov) = ion_plot_prediction_coverage(coverageNeeded);
#plt.savefig(ontargetdir + "/predicted_coverage.png"); #REMEMBER check this why

coverageNeeded = array(coverageNeeded);
meanCoveragePerBase = mean(coverageNeeded);
numBaseGT02Mean = sum(coverageNeeded>=0.2*meanCoveragePerBase);
percBaseGT02Mean = 100*(float(numBaseGT02Mean)/float(numBasesInTarget));
numBaseGT01Mean = sum(coverageNeeded>=0.1*meanCoveragePerBase);
percBaseGT01Mean = 100*(float(numBaseGT01Mean)/float(numBasesInTarget));
numBaseGT001Mean = sum(coverageNeeded>=0.01*meanCoveragePerBase);
percBaseGT001Mean = 100*(float(numBaseGT001Mean)/float(numBasesInTarget));

#meanCoveragePerBaseByChr = {}
#meanCoveragePerBaseByChr = defaultdict(float)
#numBaseGT02MeanByChr = {}
#numBaseGT02MeanByChr = defaultdict(float)
#percBaseGT02MeanByChr = {}
#percBaseGT02MeanByChr = defaultdict(float)


#coverageNeededByChr = array(coverageNeededByChr);
#meanCoveragePerBaseByChr = mean(coverageNeededByChr);
#numBaseGT02MeanByChr = sum(coverageNeededByChr>0.2*meanCoveragePerBaseByChr);
#percBaseGT02MeanByChr = 100*(float(numBaseGT02MeanByChr)/float(numBasesInTargetByChr));
#for k in coverageNeededByChr.keys():
#	coverageNeededByChr[k] = array(coverageNeededByChr[k])
#	meanCoveragePerBaseByChr[k] = mean(coverageNeededByChr[k])
#	numBaseGT02MeanByChr[k] = sum(coverageNeededByChr[k]>0.2*meanCoveragePerBaseByChr[k])
#	percBaseGT02MeanByChr[k] = 100*(float(numBaseGT02MeanByChr[k]/float(numBasesInTargetByChr[k])))

	#print 'coverageNeededByChr[k] is           ' + `k` + `coverageNeededByChr[k]` +  '    ' + `coverageNeeded`;
	#print 'meanCoveragePerBaseByChr    ' + `k` + `meanCoveragePerBaseByChr[k]` +  '    ' + `meanCoveragePerBase`;
	#print 'numBaseGT02MeanByChr[k]      ' + `k` + `numBaseGT02MeanByChr[k]` +  '    ' + `numBaseGT02Mean`;
	#print 'percBaseGT02MeanByChr[k]            ' + `k` + `percBaseGT02MeanByChr[k]` +  '    ' + `percBaseGT02Mean`; 
#	print 'meanCoveragePerBaseByChr  ' + `k` + ' ' + `meanCoveragePerBaseByChr[k]`
#	print 'numBaseGT02MeanByChr      ' + `k` + ' ' + `numBaseGT02MeanByChr[k]`
#	print 'percBaseGT02MeanByChr     ' + `k` + ' ' + `percBaseGT02MeanByChr[k]`

#99 percentile 1x and 98 percentile 10x, using the method MA showed with his Excel spreadsheet
#numNotCoveredBase = int(numBasesInTarget - coverageLineCt);
#allBaseCoverage = [];
#allBaseCoverage = [0] * numNotCoveredBase;
#depthSorted = sorted(depthCopy);
#for i in depthSorted:
#        allBaseCoverage.append(i);

#coverageNeeded9901 = estCov[0];
#coverageNeeded9820 = estCov[1];
#coverageNeeded95350 = estCov[2];

#print "time mark 11";
#print strftime("%a, %d %b %Y %X +0000", gmtime());


avgReadAcc = -1.;
#posErrSumFile = ontargetdir + "/position_error_summary.txt";
if (os.path.isfile(posErrSumFile)):
	#f = open(ontargetdir + "/position_error_summary.txt", 'r');
	f = open(posErrSumFile, 'r');
	lines = f.readlines();
	f.close();
	avgReadAcc=float(lines[5].strip().split("=")[1].strip());

if (os.path.isfile(readLenSummary)):
	f = open(readLenSummary, 'r');
	lines = f.readlines();
	f.close();
	readLenCount1=lines[0].strip().split("\t")[1].strip();
	readLenCount2=lines[1].strip().split("\t")[1].strip();
	readLenCount3=lines[2].strip().split("\t")[1].strip();
	readLenCutoff1=lines[0].strip().split("\t")[2].strip();
	readLenCutoff2=lines[1].strip().split("\t")[2].strip();
	readLenCutoff3=lines[2].strip().split("\t")[2].strip();
else:
	readLenCount1="-1";
	readLenCount2="-1";
	readLenCount3="-1";
	readLenCutoff1="-1";
	readLenCutoff2="-1";
	readLenCutoff3="-1";


#print "time mark 12";
#print strftime("%a, %d %b %Y %X +0000", gmtime());

#20130415 new requested statistics
#biased.base.file
biasedBaseNum = -1;
biasedBasePerc = -1.;
if (os.path.isfile(biasedBaseFile)):
	f = open(biasedBaseFile, 'r');
	lines = f.readlines();
	f.close();
	for line in lines:
		if 'bases with biased coverage' in line:
			biasedBaseNum = int(line.strip().split(" ")[0]);
			#print "biased base number is "
			#print biasedBaseNum
			biasedBasePerc = 100*(float(biasedBaseNum) / float(numBasesInTarget));

#read.info.file
lowMapQualReadNum = -1;
lowMapQualReadPerc = -1.;
if (os.path.isfile(readInfoFile)):
	f = open(readInfoFile, 'r');
	lines = f.readlines();
	f.close();
	lowMapQualReadNum = int(lines[-1]);
	lowMapQualReadPerc = 100*(float(lowMapQualReadNum) / float(totalReads));
#print "low map quality base num"
#print lowMapQualReadNum
#print lowMapQualReadPerc

#print "time mark 13";
#print strftime("%a, %d %b %Y %X +0000", gmtime());


fhJ.write("{\n");

str = ("Percent greater than 0.01 mean reads per base\t%2.2f%%\n" % percBaseGT001Mean);
write_to_both_files(str);
str = ("Percent greater than 0.1 mean reads per base\t%2.2f%%\n" % percBaseGT01Mean);
write_to_both_files(str);
str = ("Percent greater than 0.2 mean reads per base\t%2.2f%%\n" % percBaseGT02Mean);
write_to_both_files(str);

str= ("Percent no strand bias of all bases\t%2.2f%%\n" % percNoBaseBias);
write_to_both_files(str);
str= ("Percent no strand bias of all Amps 40~60%%\t%2.2f%%\n" % (100.*float(numNoStrandBias4060)/float(numTargets)));
write_to_both_files(str);
str= ("Percent no strand bias of all Amps 30-70%%\t%2.2f%%\n" % (100.*float(numNoStrandBias)/float(numTargets)));
write_to_both_files(str);
if (len(ampliconWithLessThan100Reads) != numTargets):
	str = ("Percent no strand bias of all Amps 30-70%% >100 reads\t%2.2f%%\n" % (100.*float(numNoStrandBias - numNoStrandBiasLess100)/(float(numTargets)-len(ampliconWithLessThan100Reads))));
else:
	str = ("Percent no strand bias of all Amps 30-70% >100 reads\t-1.00%\n");
write_to_both_files(str);
#str = ("Percent end to end read of on target reads\t%2.2f%%\n" % percEndRead);
#write_to_both_files(str);
str = ("Per base accuracy\t%2.2f%%\n" % avgReadAcc);
write_to_both_files(str);
str = ("Percent of total reads mapped to hg19\t%2.2f%%\n" % (100.*float(mappedReads)/float(totalReads)));
write_to_both_files(str);
#str = ("Percent wells with read\t%2.2f%%\n" % wellPerc);
str = ("Percent wells with read\t%2.2f%%\n" % percAddWell);
write_to_both_files(str);

str = ("Number of total reads\t%s\n" % locale.format("%d", totalReads, grouping=False));
write_to_both_files(str);
str = ("Number of mapped reads\t%s\n" % locale.format("%d", mappedReads, grouping=False));
write_to_both_files(str);
str = ("Number of targets\t%s\n" % locale.format("%d", len(targetRecords), grouping=False));
write_to_both_files(str);
str = ("Number of reads on target\t%s\n" % locale.format("%d", targetReads, grouping=False));
write_to_both_files(str);
str = ("Percent all reads on target\t%2.2f%%\n" % percAllOnTarget);
write_to_both_files(str);
str = ("Percent mapped reads on target\t%2.2f%%\n" % percMappedOnTarget);
write_to_both_files(str);
str = ("Percent reads off target\t%2.2f%%\n" % percOffTarget);
write_to_both_files(str);
str = ("Percent reads unmapped\t%2.2f%%\n" % percUnmapped);
write_to_both_files(str);
str = ("Bases in targeted reference\t%s\n" % locale.format("%d", numBasesInTarget, grouping=False));
write_to_both_files(str);
str = ("Bases covered (at least 1x)\t%s\n" % locale.format("%d", basesCovered, grouping=False));
write_to_both_files(str);
str = ("Total base reads on target\t%s\n" % locale.format("%d", totalbasereads, grouping=False));
write_to_both_files(str);
str = ("Average base coverage depth\t%s\n" % locale.format("%2.2f", avgBaseCovDepth, grouping=False));
write_to_both_files(str);
str = ("Maximum base read depth\t%s\n" % locale.format("%d", maxDepth, grouping=False));
write_to_both_files(str);
str = ("Average base read depth\t%s\n" % locale.format("%2.2f", avgBaseReadDepth, grouping=False));
write_to_both_files(str);
str = ("Std.Dev base read depth\t%s\n" % locale.format("%2.2f", sdBaseReadDepth, grouping=False));
write_to_both_files(str);

#fhC.write ("Coverage\t1x\t10x\t20x\t100x\t200x\t500x\n");
#fhC.write ("Raw target coverage\t%2.2f%%\t%2.2f%%\t%2.2f%%\t%2.2f%%\t%2.2f%%\t%2.2f%%\n" %(atLeast1x,atLeast10x,atLeast20x,atLeast100x,atLeast200x,atLeast500x));
#fhC.write ("Target coverage normalized as 100x\t%2.2f%%\t%2.2f%%\t%2.2f%%\t%2.2f%%\t%2.2f%%\t%2.2f%%\n" %(atLeast1xNorm100,atLeast10xNorm100,atLeast20xNorm100,atLeast100xNorm100,atLeast200xNorm100,atLeast500xNorm100));

str = ("Target coverage at 1x\t%2.2f%%\n" % atLeast1x);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 10x\t%2.2f%%\n" % atLeast10x);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 20x\t%2.2f%%\n" % atLeast20x);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 100x\t%2.2f%%\n" % atLeast100x);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 200x\t%2.2f%%\n" % atLeast200x);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 500x\t%2.2f%%\n" % atLeast500x);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 1x - norm 100\t%2.2f%%\n" % atLeast1xNorm100);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 10x - norm 100\t%2.2f%%\n" % atLeast10xNorm100);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 20x - norm 100\t%2.2f%%\n" % atLeast20xNorm100);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 100x - norm 100\t%2.2f%%\n" % atLeast100xNorm100);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 200x - norm 100\t%2.2f%%\n" % atLeast200xNorm100);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));
str = ("Target coverage at 500x - norm 100\t%2.2f%%\n" % atLeast500xNorm100);
write_to_both_files(str);
#fhJ.write (conv_summary_row_to_json(str));

str = ("Percent end to end read of on target reads\t%2.2f%%\n" % percEndRead);
write_to_both_files(str);
str = ("Percent forward end to end read of on target reads\t%2.2f%%\n" % percEndReadfw);
write_to_both_files(str);
str = ("Percent reverse end to end read of on target reads\t%2.2f%%\n" % percEndReadrc);
write_to_both_files(str);

str = ("Number HSM Loci >500 Reads of Priority 1 design\t%2.2f%%\n" % percHSMPri1);
write_to_both_files(str);
str = ("Number HSM Loci >500 Reads of Priority 9 design\t%2.2f%%\n" % percHSMPri9);
write_to_both_files(str);

str = ("Number of read with %s\t%s\n" % (readLenCutoff1, readLenCount1));
write_to_both_files(str);
str = ("Number of read with %s\t%s\n" % (readLenCutoff2, readLenCount2));
write_to_both_files(str);
str = ("Number of read with %s\t%s\n" % (readLenCutoff3, readLenCount3));
write_to_both_files(str);

str = ("Conversion rate at 0.2*Mean\t%2.2f%%\n" % conversionRate);
write_to_both_files(str);
str = ("Conversion rate at 0.01*Mean\t%2.2f%%\n" % conversionRate001);
write_to_both_files(str);

str = ("Base number with biased coverage (both strands >20x and one strand <3x)\t%s\n" % locale.format("%d", biasedBaseNum, grouping=False));
write_to_both_files(str);
str = ("Base percentage with biased coverage (both strands >20x and one strand <3x)\t%2.2f%%\n" % biasedBasePerc);
write_to_both_files(str);

str = ("Read number with <4 map quality\t%s\n" % locale.format("%d", lowMapQualReadNum, grouping=False));
write_to_both_files(str);
str = ("Read percentage with <4 map quality\t%2.2f%%\n" % lowMapQualReadPerc);
write_to_both_files(str);

str = ("Number of read starting >10bp from start\t%s\n" % locale.format("%d", numIn10, grouping=False));
write_to_both_files(str);
str = ("Percent of read starting >10bp from start\t%2.2f%%\n" % percReadIn10);
write_to_both_files(str);

str = ("Number of read starting ok (within 5bp) but ending early\t%s\n" % locale.format("%d", numStartOkEndEarly, grouping=False));
write_to_both_files(str);
str = ("Number of read starting ok (within 5bp) but ending early\t%2.2f%%\n" % percStartOkEndEarly);
write_to_both_files(str);


fhJ.write ("\t\"mean_covdepth\":\"%s\",\n" % locale.format("%2.2f", mu, grouping=False));
fhJ.write ("\t\"avg_num_reads_per_amp\":\"%d\",\n" % meanNumReads);
fhJ.write ("\t\"num_total_wells\":\"%d\",\n" % wellNum);
fhJ.write ("\t\"perc_wells_with_read\":\"%2.2f%%\",\n" % wellPerc);

fhJ.write ("\t\"ChrX_mean_Coverage_Per_Base_By_Chr\":\"%2.2f%%\",\n" % ChrX_meanCoveragePerBaseByChr);
fhJ.write ("\t\"ChrX_num_Base_GT_02_Mean_By_Chr\":\"%d\",\n" % ChrX_numBaseGT02MeanByChr);
fhJ.write ("\t\"ChrX_perc_Base_GT_02_Mean_By_Chr\":\"%2.2f%%\",\n" % ChrX_percBaseGT02MeanByChr);

fhJ.write ("\t\"ChrY_mean_Coverage_Per_Base_By_Chr\":\"%2.2f%%\",\n" % ChrY_meanCoveragePerBaseByChr);
fhJ.write ("\t\"ChrY_num_Base_GT_02_Mean_By_Chr\":\"%d\",\n" % ChrY_numBaseGT02MeanByChr);
fhJ.write ("\t\"ChrY_perc_Base_GT_02_Mean_By_Chr\":\"%2.2f%%\",\n" % ChrY_percBaseGT02MeanByChr);

fhJ.write ("\t\"ChrAuto_mean_Coverage_Per_Base_By_Chr\":\"%2.2f%%\",\n" % ChrAuto_meanCoveragePerBaseByChr);
fhJ.write ("\t\"ChrAuto_num_Base_GT_02_Mean_By_Chr\":\"%d\",\n" % ChrAuto_numBaseGT02MeanByChr);
fhJ.write ("\t\"ChrAuto_perc_Base_GT_02_Mean_By_Chr\":\"%2.2f%%\",\n" % ChrAuto_percBaseGT02MeanByChr);

fhJ.write ("\t\"bedfile\":\"%s\"\n" % targetname);
fhJ.write ("}\n");

fhJ.close();
fhT.close();
#fhC.close();
