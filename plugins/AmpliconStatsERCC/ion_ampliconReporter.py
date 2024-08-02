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
    print("usage: ion_ampliconReporter.py <bam.file> <bed.file> on_or_off <amp.sum.csv> <bedString> <amplicon.seq>");
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
#bedFile = sys.argv[2];
bedFileString = sys.argv[2];
checkOffOn = sys.argv[3];
ampsum = sys.argv[4];
bedString = sys.argv[5];
fhampsum = open(ampsum, 'w');

now = datetime.datetime.now()
print "Beginning bam/sam file processing: "
print now.strftime("%Y-%m-%d %H:%M:%S")

amplicons = [];
field = {};
totalNumReads = 0;
#for i in range(666):
for fn in os.listdir("."):
    if (fn.endswith(".bed") and fn.startswith("split.")):
	bf = fn;
	#print "working on bedfile %s" % bf;
	mergedStart = 0;
	mergedEnd = 0;
	mergeChr = '';

	f = open(bf, 'r');
	lines = f.readlines();
	mergedChr = lines[0].strip().split("\t")[0];
	mergedStart = lines[0].strip().split("\t")[1];
	mergedEnd = lines[-1].strip().split("\t")[2];
	
	#bedRecords = read_bed_records(bedFile);
	bedRecords = read_bed_records(bf);
	coorAmp = {};

	chromosome = defaultdict(list);


        #bed file need to be sorted
	sStarts = defaultdict(list); #sorted starts of each chrom
	sStops = defaultdict(list); #sorted stops of each chrom
	sAmplicons = defaultdict(list); #sorted amplicons of each chrom

	# setup the amplicons, 1 for each element in the bed file
	for r in range(len(bedRecords)):
		b = bedRecords[r];
		rec = b #bedRecords[b];
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

	#print "begin samtools view";
	#print now.strftime("%Y-%m-%d %H:%M:%S")
	#print " %d" % now.microsecond

	cmd = "samtools view -F 0x4 " + bamFile + " " + str(mergedChr) + ":" + str(mergedStart) + "-" + str(mergedEnd) + " | cut -f2-99";
	pipe = Popen(cmd, shell=True, stdout=PIPE).stdout

	md_patt = re.compile('\^[ACGTN]+')
	digit_patt = re.compile('\d');

	lineCt = 0;

	preChrom = '';
	preStart = 0;
	preStop = 0;
	preBestAmp = [];
	preMinDist = 0;

	
	for l in pipe:
		lineCt += 1;
		#if (lineCt % 100000 == 0):
			#now = datetime.datetime.now()
			#print "processed line %d: " % lineCt
			#print now.strftime("%Y-%m-%d %H:%M:%S")
			#print " %d" % now.microsecond

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

		(mappedLength, counts, depth) = calcMappedLengthPlus(strand, readStart, cigar);
		counts["MM"] = len(digit_patt.sub("", md_patt.sub("", md)));


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

		else:
			#print "bisect_left sort/searching %d" % len(sStarts[chrom]);
			index1 = bisect.bisect_left(sStarts[chrom], readStart);
			index11 = index1 - 1;
			index11 = max(0, index11);
			index1 = index11;

			for i in range(index1, len(sStops[chrom])):
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

	pipe.close();

now = datetime.datetime.now()

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

