#!/usr/bin/env python
import sys;
import fileinput;
import os;
from subprocess import *;

# a simple class to store suspicious regions
class TargetRegion(object):
	def __init__(self, chrom, start, stop):
		self.chrom = chrom;
		self.start = start;
		self.stop = stop;
	def __repr__(self):
		return "%s\t%d\t%d" % (self.chrom, self.start, self.stop)

# boiler plate argument processing
if len(sys.argv) != 2:
	print("usage: ion_buildOffTargetBed.py <offtarget.bam.file>");
	exit(1);

offtargetbam = sys.argv[1];

#Check files exist
if not os.path.exists(offtargetbam):
	print("BAM file does not exist");
	exit(1);


# define some constants
MIN_DEPTH = 100;
MIN_SIZE = 10;

# for every mapped on target read, set the mask = 1
cmd = "samtools depth " + offtargetbam + " | awk -v MIN_DEPTH="+ str(MIN_DEPTH) +" '{if($3>MIN_DEPTH){print}}'"

# keep an array of suspicious regions
targets = [];

# initialize the target
t = TargetRegion("none", -1, -1);

#open the pipe and process each line
pipe = Popen(cmd, shell=True, stdout=PIPE).stdout
for l in pipe:
	l = l.strip().split();
	chrom = l[0];
	pos = int(l[1]);
	if not chrom==t.chrom:
		if (t.stop-t.start) >= MIN_SIZE:
			targets.append(t);
		t = TargetRegion(chrom, pos, pos);	
	else:
		if pos == (t.stop+1):
			t.stop = pos;
		else:
			if (t.stop-t.start) >= MIN_SIZE:
				targets.append(t);
			t = TargetRegion(chrom, pos, pos);
pipe.close();

# print out the targets in a bed format
if not t.chrom == 'none':
	targets.append(t);
for t in targets:
	print(t);
