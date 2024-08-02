#!/usr/bin/env python

# A small script to output the summary of reads on a per-chromosome basis.
# Andrew Hanna, 2012, IonTorrent.

import sys;
import fileinput;
from collections import defaultdict;
from subprocess import *;
from ampliconUtils import *;
from numpy import *;
import datetime;

# A container class for chromosome numbers
class ChromStats(object):
	def __init__(self, name):
		self.name = name;
		self.numFwd = 0;
		self.numRev = 0;
		self.depth = [];
	def __repr__(self):
		return self.name + "\t" + str(self.numFwd) + "\t" + str(self.numRev);

# sort the list l in natural order
def natural_sort(l): 
	convert = lambda text: int(text) if text.isdigit() else text.lower() 
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
	return sorted(l, key = alphanum_key)

# parse the input parameters

if len(sys.argv) != 4:
    print("usage: ion_chromosomeReporter.py <bam.file> <bed.file> <coverage.file>");
    exit(0);
bamFile = sys.argv[1];
bedFile = sys.argv[2];
coverageFile = sys.argv[3];

# 1. parse the bam file with reads aligned to the genome
# 2. only consider mapped & on-target reads
# 3. maintain a count of the number of reads per chromosome.

now = datetime.datetime.now()
print >> sys.stderr, "Begin to samtool view: "
print >> sys.stderr, now.strftime("%Y-%m-%d %H:%M:%S")
cmd = "samtools view -F 0x4 -L" + bedFile + " " + bamFile + " | cut -f2-3";
chromosomes = defaultdict(lambda: ChromStats("Unknown"));
pipe = Popen(cmd, shell=True, stdout=PIPE).stdout
for l in pipe:
	l = l.strip();
	tokens = l.split();
	strand = int(tokens[0]);
	chrom = tokens[1];
	c = chromosomes[chrom];
	c.name = chrom;
	if strand == 0:
		c.numFwd += 1;
	elif strand == 16:
		c.numRev += 1;
	chromosomes[chrom] = c;
pipe.close();
now = datetime.datetime.now()
print >> sys.stderr, "finished samtool view: "
print >> sys.stderr, now.strftime("%Y-%m-%d %H:%M:%S")


# 4. From the coverage file, get the depth counts for each base covered over the chromosome
for line in fileinput.input(coverageFile):
	line = line.strip();
	tokens = line.split();
	c = chromosomes[tokens[0]];
	c.name = tokens[0];
	c.depth.append(int(tokens[2]));

# 5. print out the table to standard out.

print "Chrom,BasesCovered,TotalBases,FwdReads,RevReads,TotalReads,1x,10x,100x,500x,1000x";
for chrom in natural_sort(chromosomes.keys()):
	c = chromosomes[chrom]; 
	d = array(c.depth);
	cov1x = sum(d>=1);
	cov10x = sum(d>=10);
	cov100x = sum(d>=100);
	cov500x = sum(d>=500);
	cov1000x = sum(d>=1000);
	print "%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d" % (c.name, len(c.depth), sum(c.depth), c.numFwd, c.numRev, c.numFwd + c.numRev, cov1x, cov10x, cov100x, cov500x, cov1000x);
