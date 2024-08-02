import fileinput;
import sys;
from collections import defaultdict;
from math import *;
import re;
import operator;

class Amplicon(object):
	def __init__(self, name, numfwd, numrev, seq, numfwde2e, numreve2e):
	#def __init__(self, name, numfwd, numrev, seq, numfwde2e, numreve2e, fwduprimer, revuprimer):
		self.name = name;
		self.numfwd = numfwd;
		self.numrev = numrev;
		self.seq = seq;
		self.total = numfwd + numrev;
		self.numfwde2e = numfwde2e;
		self.numreve2e = numreve2e;
		self.totale2e = numfwde2e + numreve2e;
		self.counts = defaultdict(int);
		self.depth = defaultdict(int);
#		self.depthF = defaultdict(int);
#		self.depthR = defaultdict(int);
                #starts > 10bp into bed region, i.e., starts late
		self.numfwdin10 = 0;
		self.numrevin10 = 0;
		self.numfwdstartokayendearly = 0; #starts <= 5bp into bed region starting point, end >5bp earlier
		self.numrevstartokayendearly = 0; 

		self.chrom = '';
		self.insStart = 0;
		self.insStop = 0;

		self.insIQR = 0.;
		self.insAVG = 0.;
#		self.insIQRavgF = 0.;
#		self.insIQRavgR = 0.;

		self.insAvg02Perc = 0.;

#		self.fwduprimer = fwduprimer;
#		self.revuprimer = revuprimer;
		self.fwduprimer = '';
		self.revuprimer = '';
		self.fwdprimer = '';
		self.revprimer = '';
		self.note = '';
		self.field = {}

	def setChrom(self, chrom):
		self.chrom = chrom;
	def setInsStart(self, insStart):
		self.insStart = insStart;
	def setInsStop(self, insStop):
		self.insStop = insStop;

	def addCounts(self, cnts):
		for c in cnts:
			self.counts[c] += cnts[c];
	def getTotal(self):
		return self.numfwd + self.numrev;

	def getTotalE2E(self):
		return self.numfwde2e + self.numreve2e;

	def getGC(self):
		if len(self.seq) == 0:
			return 0;
		GC_count = float(self.seq.count('G') + self.seq.count('C'));
		return GC_count/float(len(self.seq));

	def getBias(self):
		total = float(self.getTotal());
		if total == 0:
			return 0;
		else:
			return float(self.numfwd)/total;

#	def setFwdUPrimer(self, fwduprimer):
#		self.fwduprimer = fwduprimer;
#	def setRevUPrimer(self, revuprimer):
#		self.revuprimer = revuprimer;
#	def getFwdUPrimer(self):
#		return self.fwduprimer;
#	def getRevUPrimer(self):
#		return self.revuprimer;


	def __repr__(self):
		DELIM = "\t";
		return self.name + DELIM + self.seq + DELIM + str(self.numfwd) + DELIM + str(self.numrev) + DELIM + str(self.getTotal()) + DELIM + str(self.getBias()) ;

class FastaRecord(object):
	def __init__(self, title, sequence):
		self.title = title
		self.seq = sequence
	def __repr__(self):
		return ">" + self.title + "\n" + self.seq;

class TargetRegion(object):
#	def __init__(self, title, start, stop, chrom):
	def __init__(self, title, start, stop, chrom, note):
		self.title = title;
		self.chrom = chrom;
		self.chrStart = start;
		self.chrStop = stop;
                self.numfwd = 0;
                self.numrev = 0;
                self.depth = [];
		self.note = note;
	def contains(self, chrom, readStart, readStop):
		if (self.chrom==chrom):
			if (readStop<self.chrStart) | (readStart>self.chrStop):
				return False;
			return True;
#				return True;
		return False;

	def __repr__(self):
#		DELIM = "\t";
		DELIM = ",";
		return self.chrom + DELIM + str(self.chrStart) + DELIM + str(self.chrStop) + DELIM + self.title;
#                return self.title + DELIM + self.chrom + DELIM + str(self.chrStart) + DELIM + str(self.chrStop);
        def setNumFwd(self, fwd):
                self.numfwd = fwd;
        def setNumRev(self, rev):
                self.numrev = rev;
        def getLength(self):
                return self.chrStop - self.chrStart;

def calcAccuracy(counts):
	total_error = float(counts['D'] + counts['I'] + counts['MM']);
	total_base = float(counts['M'] + counts['I']);
	if (total_base == 0):
		acc = -1.;
	else:
		acc = 100.*(1. - (total_error/total_base));
	return acc;
	#return(acc, total_error, total_base);


def calcMappedLength(cigar):
	counts = defaultdict(int);
	pcigar = defaultdict(int);
	pattern = re.compile('([0-9]*)([DMIX=])');
	total = 0;	
	for n, c in pattern.findall(cigar):
		counts[c] += int(n);
		if not c=="I":
			total += int(n);
	return (total, counts);

def calcMappedLengthPlus(strand, rstart, cigar):
	counts = defaultdict(int);
	pcigar = defaultdict(int);
	depth = defaultdict(int);
#	depthF = defaultdict(int);
#	depthR = defaultdict(int);

	pattern = re.compile('([0-9]*)([DMIX=])');
	total = 0;
	mark = rstart;
	for n, c in pattern.findall(cigar):
		if (c == 'M' or c == 'X' or c == '='):
			for i in range(int(n)):
				depth[mark] = 1;
#				if strand == 0:
#					depthF[mark] = 1;
#				elif strand == 16:
#					depthR[mark] = 1;
				mark += 1;
		elif (c == 'D'):
			for i in range(int(n)):
				depth[mark] = 0;
#				if strand == 0:
#					depthF[mark] = 0;
#				elif strand == 16:
#					depthR[mark] = 0;
				mark += 1;
#		elif (c == 'I'):
			#do nothing
			
		counts[c] += int(n);
		if not c=="I":
			total += int(n);
#	return (total, counts, depth, depthF, depthR);
	return (total, counts, depth);



def read_bed_records(fname):
	records = defaultdict(TargetRegion)
	inRec = False;
	seq = "";
	title = "";
	num = 1;
	note = '';
	ct = 0;
	for line in fileinput.input(fname):
		if not line.startswith('track'):
			line = line.strip();
			tokens = line.split();
			chrom = tokens[0];
			start = int(tokens[1]);
			stop = int(tokens[2]);
			if (len(tokens)>3):
				title = tokens[3];
				if (len(tokens)>6):
					note = tokens[6]
			else:
				title = "unknown_" + str(num);
				num +=1;
			#records[title] = TargetRegion(title, start, stop, chrom);
			#records[title] = TargetRegion(title, start, stop, chrom, note);
			records[ct] = TargetRegion(title, start, stop, chrom, note);
			ct += 1;
	return records

def read_chr_bed_records(fname, chrName):
	records = defaultdict(TargetRegion)
	inRec = False;
	seq = "";
	title = "";
	num = 1;	
	for line in fileinput.input(fname):
		if not line.startswith('track'):
			line = line.strip();
			tokens = line.split();
			chrom = tokens[0];

			if chrom.lower() == chrName.lower():
				start = int(tokens[1]);
				stop = int(tokens[2]);
				if (len(tokens)>3):
					title = tokens[3];
				else:
					title = "unknown_" + str(num);
					num +=1;
				records[title] = TargetRegion(title, start, stop, chrom);
	return records

def read_fasta_records(fname):
	records = defaultdict(FastaRecord);
	inRec = False;
	seq = "";
	title = "";
	ct = 0;
	for line in fileinput.input(fname):
		if line.startswith('>'):
			if inRec:
				records[title] = FastaRecord(title, seq);
				#records[ct] = FastaRecord(title, seq);
			inRec = True;
			seq = "";
			title = line[1:].rstrip();
			ct += 1;
		else:
			line = line.upper();
			seq += line.rstrip();
	if inRec:
		records[title] = FastaRecord(title, seq);
		#records[ct] = FastaRecord(title, seq);
	return records




#def bisect_left(a, x, lo=0, hi=None):
	"""Return the index where to insert item x in list a, assuming a is sorted.
     
        The return value i is such that all e in a[:i] have e < x, and all e in
        a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
        insert just before the leftmost x already there.
    
        Optional args lo (default 0) and hi (default len(a)) bound the
        slice of a to be searched.
        """
   
#        if lo < 0:
#            raise ValueError('lo must be non-negative')
#        if hi is None:
#            hi = len(a)
#        while lo < hi:
#            mid = (lo+hi)//2
#            if a[mid] < x: lo = mid+1
#            else: hi = mid
#        return lo

