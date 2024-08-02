#!/usr/bin/env python
import sys;
import os;
import fileinput;
from collections import defaultdict;
from subprocess import *;
from ampliconUtils import *;
import random as rnd;
from numpy import *;
from scipy.stats import *;
import datetime;
import bisect;

def isCloseToStartOfAmp(ampPos, readStart, direc):
	if direc==0:
		if readStart<(ampPos+5):
			return True;
	elif direc==16:
		if readStart>(ampPos-5):
			return True;
	return False;



if len(sys.argv) < 6:
    print("usage: ion_ampliconReporter.py <bam.file> <bed.file> on_or_off <mid.F0x4.sam> <amp.sum.csv> <amplicon.fasta>");
    exit(0);
if len(sys.argv) == 7:
	ampliconFile = sys.argv[6];
else:
	ampliconFile = "";

if ((len(ampliconFile)>0) & (os.path.isfile(ampliconFile))):
	ampliconSequences = read_fasta_records(ampliconFile);
else:
	ampliconSequences = defaultdict(lambda: FastaRecord("Unknown", ""));


bamFile = sys.argv[1];
bedFile = sys.argv[2];
checkOffOn = sys.argv[3];
midF0x4Sam = sys.argv[4];
ampsum = sys.argv[5];
fhampsum = open(ampsum, 'w');

now = datetime.datetime.now()
print "Beginning bed records processing: "
print now.strftime("%Y-%m-%d %H:%M:%S")

bedRecords = read_bed_records(bedFile);
coorAmp = {};

chromosome = defaultdict(list);
field = {};
amplicons = [];

#bed file need to be sorted
sStarts = defaultdict(list); #sorted starts of each chrom
sStops = defaultdict(list); #sorted stops of each chrom
sAmplicons = defaultdict(list); #sorted amplicons of each chrom

# setup the amplicons, 1 for each element in the bed file
for r in range(len(bedRecords)):
	#print r
#for b in bedRecords:
	b = bedRecords[r];
	#print b
	rec = b #bedRecords[b];
	#if (rec.chrom == "chrY"):
	#	print "r is %d" % r
	#	print b
	seq = ampliconSequences[b.title];
	amplicon = Amplicon(b.title, 0, 0, seq.seq, 0, 0);
	amplicon.target = rec;
	amplicon.numFiltFwd = 0;
	amplicon.numFiltRev = 0;
	amplicon.numSnapFwd = 0;
	amplicon.numSnapRev = 0;
	chromosome[rec.chrom].append(amplicon);
	amplicons.append(amplicon);

        amplicon.chrom = rec.chrom;
        amplicon.insStart = rec.chrStart;
        amplicon.insStop = rec.chrStop;

	sStarts[rec.chrom].append(rec.chrStart);
	sStops[rec.chrom].append(rec.chrStop);
	sAmplicons[rec.chrom].append(amplicon);
	#coor = "%d_%d" % (rec.chrStart, rec.chrStop)
	coor = "%s_%d_%d" % (rec.chrom, rec.chrStart, rec.chrStop)
	coorAmp[coor] = amplicon

	if not ';' in rec.note:
		continue;
	if not 'on' in checkOffOn:
		continue;

	amplicon.note = rec.note;
	tokens = rec.note.split(';');
	for kv in tokens:
		fs = kv.split('=');
		try:
			field[fs[0]] = fs[1];
			(amplicon.field)[fs[0]] = fs[1]
		except:
			field[fs[0]] = 'NA';
			(amplicon.field)[fs[0]] = 'NA'

		if (fs[0] == 'FWD_UPRIMER'):
			try:
				amplicon.fwduprimer = fs[1];
				amplicon.fwdprimer = fs[1].replace('U', 'T');
				amplicon.fwdprimer = fs[1].replace('u', 't');
			except:
				amplicon.fwduprimer = 'NA';
		if (fs[0] == 'REV_UPRIMER'):
			try:
				amplicon.revuprimer = fs[1];
				amplicon.revprimer = fs[1].replace('U', 'T');
				amplicon.revprimer = fs[1].replace('u', 't');
			except:
				amplicon.revuprimer = 'NA';


#print sStarts["chrY"]
#print sStops["chrY"]
#print sAmplicons["chrY"]
#for c in sStarts:
#	sStarts[c] = sort(sStarts[c]);
#print sStarts["chr1"]

now = datetime.datetime.now()
print "Beginning bam/sam file processing: "
print now.strftime("%Y-%m-%d %H:%M:%S")


cmd = "samtools view -F 0x4 -L" + bedFile + " " + bamFile + " | cut -f2-99";
#cmd = "samtools view -F 0x4 -L" + bedFile + " -S " + midF0x4Sam + " | cut -f2-99";
pipe = Popen(cmd, shell=True, stdout=PIPE).stdout

md_patt = re.compile('\^[ACGTN]+')
digit_patt = re.compile('\d');

totalNumReads = 0;
lineCt = 0;

preChrom = '';
preStart = 0;
preStop = 0;
preBestAmp = [];
preMinDist = 0;

for l in pipe:
	lineCt += 1;
	if (lineCt % 1000000 == 0):
		now = datetime.datetime.now()
		print "processed line %d: " % lineCt
		print now.strftime("%Y-%m-%d %H:%M:%S")
		print " %d" % now.microsecond
		
	l = l.strip();
	tokens = l.split()
	#tokens.pop(0)
	strand = int(tokens[0]);
	chrom = tokens[1];
	readStart = int(tokens[2]);
	quality = int(tokens[3]);
	cigar = tokens[4];
	md = '';

        # extract the MD
        for n in range(10, len(tokens)):
            if( tokens[n].startswith('MD') ):
		    md = tokens[n];
		    md = md.split(":")[2];
		    break;

	#(mappedLength, counts) = calcMappedLength(cigar);
	#(mappedLength, counts, depth, depthF, depthR) = calcMappedLengthPlus(strand, readStart, cigar);
	(mappedLength, counts, depth) = calcMappedLengthPlus(strand, readStart, cigar);
	counts["MM"] = len(digit_patt.sub("", md_patt.sub("", md)));

	#if (lineCt % 1000 == 0):
	#	now = datetime.datetime.now()
	#	print "calcMappedLengthPlus finished for line %d: " % lineCt
	#	print now.strftime("%Y-%m-%d %H:%M:%S")
	#	print " %d" % now.microsecond


	readStop = readStart + mappedLength;
	found = False;

	minDist = 0;
	bestAmp = [];


	if (chrom == preChrom and readStart == preStart and readStop == preStop):
		preChrom = chrom;
		preStart = readStart;
		preStop = readStop;
		bestAmp = preBestAmp;
		minDist = preMinDist;
		#print "gets the same one as preivous: "
		#print bestAmp
	else:
		index1 = bisect.bisect_left(sStarts[chrom], readStart);
		index11 = index1 - 1;
		index11 = max(0, index11);
		index1 = index11;
        	#print "find bisect_left index1 %d between position %d and %d suit for %s and %d" % (index1, sStarts[chrom][index1 - 1], sStarts[chrom][index1], chrom, readStart);
		#index2 = bisect.bisect_right(sStops[chrom], readStop);
		#index22 = index2 + 1;
		#index22 = min(index22, (len(sStops[chrom])));
		#index2 = index22;

		#if (index1 >= 1):
		#	beforeIndex1 = index1 - 1;
		#	for j in range(0, beforeIndex1):
		#		if (sStarts[chrom][j] != sStarts[chrom][index1]):
		#			sStarts[chrom].pop(0);
		#			sStops[chrom].pop(0);

        	#print "find bisect_left index2 %d between position %d and %d suit for %s and %d" % (index2, sStops[chrom][index2 - 1], sStops[chrom][index2], chrom, readStop);
		#for i in range(index1, index2):
		for i in range(index1, len(sStops[chrom])):
        	#for i in range(0, (len(sStops[chrom]))):
			#coor = "%d_%d" % (sStarts[chrom][i], sStops[chrom][i]);
			coor = "%s_%d_%d" % (chrom, sStarts[chrom][i], sStops[chrom][i]);
			a = coorAmp[coor];
			if (a.target.chrStart > readStop):
				break;
			if (a.target.contains(chrom, readStart, readStop)):
				s = max(readStart, a.target.chrStart);
				e = min(readStop, a.target.chrStop);
				dist = e-s;
				if (dist >0):
					if (dist>minDist):
						minDist = dist;
						bestAmp = [a];
					elif (dist == minDist):
						bestAmp.append(a);
				else:
					break;
					

		preChrom = chrom;
		preStart = readStart;
		preStop = readStop;
		preBestAmp = bestAmp;
		preMinDist = minDist;


	if minDist > 0:
		a = bestAmp[rnd.randint(0, len(bestAmp)-1)];
		a.addCounts(counts);
		found = True;
		if strand==0:
			a.numfwd +=1;
			if readStart<=(a.target.chrStart+5):
				a.numFiltFwd += 1;
                                #if (readStop >= a.target.chrStop-5 and readStop <= a.target.chrStop+5):
                                if (readStart >= a.target.chrStart-5 and readStop-1 >= a.target.chrStop-5 and readStop-1 <= a.target.chrStop+5):
                                    a.numSnapFwd += 1;
			if (readStart <= a.target.chrStart and readStop >= a.target.chrStop):
				a.numfwde2e += 1;
			if (readStart >= a.target.chrStart + 10):
				a.numfwdin10 += 1;
			#starts <= 5bp into bed region starting point, end >5bp earlier
			if (readStart <= a.target.chrStart + 5 and readStop < a.target.chrStop - 5):
				a.numfwdstartokayendearly

		elif strand==16:
			a.numrev +=1;
			if (readStop-1)>=(a.target.chrStop-5):
				a.numFiltRev += 1;
                                if (readStop-1 <= a.target.chrStop+5 and readStart >= a.target.chrStart-5 and readStart <= a.target.chrStart+5):
                                    a.numSnapRev += 1;
			if (readStart <= a.target.chrStart and readStop >= a.target.chrStop):
				a.numreve2e += 1;
			if (readStart >= a.target.chrStart + 10):
				a.numrevin10 += 1;
			#starts <= 5bp into bed region starting point, end >5bp earlier
			if (readStart <= a.target.chrStart + 5 and readStop < a.target.chrStop - 5):
				a.numrevstartokayendearly


		totalNumReads +=1;

                for key in depth:
                    a.depth[key] += 1;
#                for key in depthF:
#                    a.depthF[key] += 1;
#                for key in depthR:
#                    a.depthR[key] += 1;

	#if (lineCt % 1000 == 0):
	#	now = datetime.datetime.now()
	#	print "bestAmp finished for line %d: " % lineCt
	#	print now.strftime("%Y-%m-%d %H:%M:%S")
	#	print " %d" % now.microsecond

pipe.close();

now = datetime.datetime.now()
print "finished all lines %d: " % lineCt
print now.strftime("%Y-%m-%d %H:%M:%S")
print " %d" % now.microsecond


for amp in amplicons:
    #amp.insIQR = -1;
    #amp.insAVG = -1;
    #amp.insAvg02Perc = -1;
    #check: comment the caclualtion to see if how much time it consumes
    #print "amplicon_name_is " + str(amp.name) + " " + str(amp.insStart) + " " + str(amp.insStop);
    insDepth = defaultdict(int);
    #insDepthF = defaultdict(int);
    #insDepthR = defaultdict(int);
    for k in range(amp.insStart, amp.insStop + 1):
        if (amp.depth[k]):
            insDepth[k] = amp.depth[k];
        else:
            insDepth[k] = 0;

    work = [];
    for k, v in insDepth.iteritems():
        work.append(v);
    #print work;

    if len(work) == 0:
        amp.insIQR = 0.;
        amp.insAVG = 0.;
    else:
        q1 = scoreatpercentile(work[:], 25);
        q3 = scoreatpercentile(work[:], 75);
        avg = mean(work[:]);
        iqr = (q3 - q1);
        amp.insIQR = iqr;
        amp.insAVG = avg;
    avg02 = 0;
    for k, v in insDepth.iteritems():
        if (v >= avg):
            avg02 += 1;
    amp.insAvg02Perc = float(avg02) * 100 / (amp.insStop - amp.insStart + 1); 
    if len(work) == 0:
        amp.insAvg02Perc = -1.;

ampReads = [];
for amp in amplicons:
	ampReads.append(amp.getTotal());
mu = mean(ampReads);


DELIM = "\t";

#print "fields"
#print field.keys()

fieldheader = '';
if (len(field.keys())>0):
	header = []
	for key in field.keys():
		fieldheader += ',%s' % key;
		header.append(key);
#	fhampsum.write ("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s\n" % ("AmpliconId","NumFwd","NumRev","Total","StrandBias","PercTotalReads","NumFwdFilt","NumRevFilt","TotalFilt","StrandBiasFilt","PercTotalReadsFilt","Pass(>0.2xMeanReads)","Num1FwdEnd2End","NumRevEnd2End","NumEnd2End","Amp GC", "Sequence", "Matches", "Inserts", "Deletions", "Mismatches", "Accuracy","NumFwdFiltSnapSize","NumRevFiltSnapSize","InsDepthIQR","InsDepthAVG","InsDepthAvg02Perc","Pass(>0.01xMeanReads)",fieldheader));
	fhampsum.write ("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s%s\n" % ("AmpliconId","NumFwd","NumRev","Total","StrandBias","PercTotalReads","NumFwdFilt","NumRevFilt","TotalFilt","StrandBiasFilt","PercTotalReadsFilt","Pass(>0.2xMeanReads)","Num1FwdEnd2End","NumRevEnd2End","NumEnd2End","Amp GC", "Sequence", "Matches", "Inserts", "Deletions", "Mismatches", "Accuracy","NumFwdFiltSnapSize","NumRevFiltSnapSize","InsDepthIQR","InsDepthAVG","InsDepthAvg02Perc","Pass(>0.01xMeanReads)","NumRead10bpIn","NumStartOkEndEarly",fieldheader));

else:
#	fhampsum.write ("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % ("AmpliconId","NumFwd","NumRev","Total","StrandBias","PercTotalReads","NumFwdFilt","NumRevFilt","TotalFilt","StrandBiasFilt","PercTotalReadsFilt","Pass(>0.2xMeanReads)","Num1FwdEnd2End","NumRevEnd2End","NumEnd2End","Amp GC", "Sequence", "Matches", "Inserts", "Deletions", "Mismatches", "Accuracy","NumFwdFiltSnapSize","NumRevFiltSnapSize","InsDepthIQR","InsDepthAVG","InsDepthAvg02Perc","Pass(>0.01xMeanReads)"));
	fhampsum.write ("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % ("AmpliconId","NumFwd","NumRev","Total","StrandBias","PercTotalReads","NumFwdFilt","NumRevFilt","TotalFilt","StrandBiasFilt","PercTotalReadsFilt","Pass(>0.2xMeanReads)","Num1FwdEnd2End","NumRevEnd2End","NumEnd2End","Amp GC", "Sequence", "Matches", "Inserts", "Deletions", "Mismatches", "Accuracy","NumFwdFiltSnapSize","NumRevFiltSnapSize","InsDepthIQR","InsDepthAVG","InsDepthAvg02Perc","Pass(>0.01xMeanReads)","NumRead10bpIn","NumStartOkEndEarly"));

#print "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % ("AmpliconId","NumFwd","NumRev","Total","StrandBias","PercTotalReads","NumFwdFilt","NumRevFilt","TotalFilt","StrandBiasFilt","PercTotalReadsFilt","Pass(>0.2xMeanReads)","Num1FwdEnd2End","NumRevEnd2End","NumEnd2End","Amp GC", "Sequence", "Matches", "Inserts", "Deletions", "Mismatches", "Accuracy","NumFwdFiltSnapSize","NumRevFiltSnapSize","InsDepthIQR","InsDepthAVG","InsDepthAvg02Perc","Pass(>0.01xMeanReads)");

ampReads = [];
for amp in amplicons:
	passFlag = 0; #original conversion rate >0.2mu
	passFlag001 = 0;

	numIn10 = amp.numfwdin10 + amp.numrevin10;
	numStartOkEndEarly = amp.numfwdstartokayendearly + amp.numrevstartokayendearly

	if amp.getTotal()>(0.2*mu):
		passFlag = 1;
	if amp.getTotal()>(0.01*mu):
		passFlag001 = 1;
	ampReads.append(amp.getTotal());
	totalFilt = amp.numFiltFwd + amp.numFiltRev;
	if totalFilt == 0:
		biasFilt = 0;
	else:
		biasFilt = float(amp.numFiltFwd)/float(totalFilt)
	if totalNumReads == 0:
		totalPerc = 0;
		totalPercFilt = 0;
	else:
		totalPerc = 100.*float(amp.getTotal())/float(totalNumReads);
		totalPercFilt = 100.*float(totalFilt)/float(totalNumReads);

	counts = amp.counts;
	

	if (len(field.keys())>0):
		fieldvalue = '' #"%s,%s" % (amp.fwdprimer, amp.revprimer);
		for k in header:
			try:
				fieldvalue += ",%s" % (amp.field)[k]
			except:
				fieldvalue += ",NA"
#		fhampsum.write ("%s,%d,%d,%d,%2.3f,%2.5f,%d,%d,%d,%2.3f,%2.5f,%d,%d,%d,%d,%2.2f,%s,%d,%d,%d,%d,%2.2f,%d,%d,%2.3f,%2.3f,%2.3f,%d%s\n" % (amp.name, amp.numfwd, amp.numrev, amp.getTotal(), amp.getBias(), totalPerc, amp.numFiltFwd, amp.numFiltRev, totalFilt, biasFilt, totalPercFilt, passFlag, amp.numfwde2e, amp.numreve2e, amp.getTotalE2E(), amp.getGC(), amp.seq, counts['M'], counts['I'], counts['D'], counts['MM'], calcAccuracy(counts), amp.numSnapFwd, amp.numSnapRev,amp.insIQR,amp.insAVG,amp.insAvg02Perc, passFlag001, fieldvalue));
		fhampsum.write ("%s,%d,%d,%d,%2.3f,%2.5f,%d,%d,%d,%2.3f,%2.5f,%d,%d,%d,%d,%2.2f,%s,%d,%d,%d,%d,%2.2f,%d,%d,%2.3f,%2.3f,%2.3f,%d,%d,%d%s\n" % (amp.name, amp.numfwd, amp.numrev, amp.getTotal(), amp.getBias(), totalPerc, amp.numFiltFwd, amp.numFiltRev, totalFilt, biasFilt, totalPercFilt, passFlag, amp.numfwde2e, amp.numreve2e, amp.getTotalE2E(), amp.getGC(), amp.seq, counts['M'], counts['I'], counts['D'], counts['MM'], calcAccuracy(counts), amp.numSnapFwd, amp.numSnapRev,amp.insIQR,amp.insAVG,amp.insAvg02Perc, passFlag001, numIn10, numStartOkEndEarly, fieldvalue));
	else:
		fhampsum.write ("%s,%d,%d,%d,%2.3f,%2.5f,%d,%d,%d,%2.3f,%2.5f,%d,%d,%d,%d,%2.2f,%s,%d,%d,%d,%d,%2.2f,%d,%d,%2.3f,%2.3f,%2.3f,%d,%d,%d\n" % (amp.name, amp.numfwd, amp.numrev, amp.getTotal(), amp.getBias(), totalPerc, amp.numFiltFwd, amp.numFiltRev, totalFilt, biasFilt, totalPercFilt, passFlag, amp.numfwde2e, amp.numreve2e, amp.getTotalE2E(), amp.getGC(), amp.seq, counts['M'], counts['I'], counts['D'], counts['MM'], calcAccuracy(counts), amp.numSnapFwd, amp.numSnapRev,amp.insIQR,amp.insAVG,amp.insAvg02Perc, passFlag001, numIn10, numStartOkEndEarly));

#	print "%s,%d,%d,%d,%2.3f,%2.5f,%d,%d,%d,%2.3f,%2.5f,%d,%d,%d,%d,%2.2f,%s,%d,%d,%d,%d,%2.2f,%d,%d,%2.3f,%2.3f,%2.3f,%d,%s,%s" % (amp.name, amp.numfwd, amp.numrev, amp.getTotal(), amp.getBias(), totalPerc, amp.numFiltFwd, amp.numFiltRev, totalFilt, biasFilt, totalPercFilt, passFlag, amp.numfwde2e, amp.numreve2e, amp.getTotalE2E(), amp.getGC(), amp.seq, counts['M'], counts['I'], counts['D'], counts['MM'], calcAccuracy(counts), amp.numSnapFwd, amp.numSnapRev,amp.insIQR,amp.insAVG,amp.insAvg02Perc, passFlag001, amp.fwduprimer, amp.revuprimer);

fhampsum.close();

