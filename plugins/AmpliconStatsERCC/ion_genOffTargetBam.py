#!/usr/bin/env python
import sys;
import fileinput;
from subprocess import *;
from numpy import *;
import os;

# boiler plate argument processing
if len(sys.argv) != 3:
	print("usage: ion_genOffTargetBam.py <bam.file> <bed.file>");
	exit(1);

bamname = sys.argv[1];
bedname = sys.argv[2];

#Check files exist
if not os.path.exists(bamname):
	print("BAM file does not exist");
	exit(1);

if not os.path.exists(bedname):
	print("BED file does not exist");
	exit(1);


# make a 2D mask array large enough for a 318 chip
#ontargetwells = zeros((100000, 100000));
ontargetwells = {}

# for every mapped on target read, set the mask = 1
cmd = "samtools view -F 0x4 -L " + bedname + " " + bamname;
pipe = Popen(cmd, shell=True, stdout=PIPE).stdout

for l in pipe:
	l = l.strip().split();
	strand = int(l[1]);
	seq = l[3];
	t = l[0].split(':');
	x = int(t[1]);
	y = int(t[2]);
	s = '';
	s += '%d' % x;
	s += '_%d' % y;
	ontargetwells[s] = 1;
	#ontargetwells[x][y] = 1;
pipe.close();


# for every mapped read, check the mask and if 0, print the read
cmd = "samtools view -F 0x4 -H " + bamname;
pipe = Popen(cmd, shell=True, stdout=PIPE).stdout
for line in pipe:
	print(line.strip());
pipe.close();

cmd = "samtools view -F 0x4 " + bamname;
pipe = Popen(cmd, shell=True, stdout=PIPE).stdout

for line in pipe:
	line = line.strip();
	l = line.split();
	strand = int(l[1]);
	seq = l[3];
	t = l[0].split(':');
	x = int(t[1]);
	y = int(t[2]);
	s = '';
	s += '%d' % x;
	s += '_%d' % y;
	if s in ontargetwells:
		a = 1;
	else:
		print(line);
	#if ontargetwells[x][y] == 0:
	#	print(line);
pipe.close();



