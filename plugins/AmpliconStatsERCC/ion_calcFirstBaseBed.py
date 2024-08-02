#!/usr/bin/env python

from subprocess import *;
import sys;
import re;

def calcMappedLength(cigar):
	pattern = re.compile('([0-9]*)([DMIX=])');
	total = 0;	
	for n, c in pattern.findall(cigar):
		if not c=="I":
			total += int(n);
	return (total);


bamname = sys.argv[1];
direction = sys.argv[2];


# for every mapped on target read, set the mask = 1
cmd = "samtools view -" + direction + "0x10 " + bamname;
pipe = Popen(cmd, shell=True, stdout=PIPE).stdout
for l in pipe:
	l = l.strip().split();
	strand = int(l[1]);
	chrom = l[2];
	startPos = int(l[3]);
	seq = l[3];
	cigar = l[5];
	if (strand == 0):
		firstBasePos = startPos;
		if (firstBasePos - 2 > 0):
			print "%s\t%d\t%d\t%s\t0\t+" % (chrom, firstBasePos-2, firstBasePos-1, "a");
	elif (strand == 16):
		firstBasePos = startPos + calcMappedLength(cigar);
		if (firstBasePos - 1 > 0):
			print "%s\t%d\t%d\t%s\t0\t-" % (chrom, firstBasePos-1, firstBasePos, "b");
pipe.close();




