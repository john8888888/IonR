#!/usr/bin/env python
from collections import defaultdict;
import fileinput;
import sys;
from matplotlib import use
use("Agg");
import matplotlib.pyplot as plt;
from numpy import *;

class TargetRegion(object):
	def __init__(self, title, start, stop, chrom):
		self.title = title;
		self.chrom = chrom;
		self.chrStart = start;
		self.chrStop = stop;
		self.depth = [];
		self.depthF = [];
		self.depthR = [];
	def __repr__(self):
		DELIM = ",";
		return self.title + DELIM + self.chrom + DELIM + str(self.chrStart) + DELIM + str(self.chrStop);



def read_bed_records(fname):
	records = defaultdict(TargetRegion)
	inRec = False;
	seq = "";
	title = "";
	for line in fileinput.input(fname):
		line = line.strip();
		tokens = line.split();
		chrom = tokens[0];
		start = float(tokens[1]);
		stop = float(tokens[2]);
		if (len(tokens)>3):
			title = tokens[3];
		else:
			title = "unknown";
		records[title] = TargetRegion(title, start, stop, chrom);
	return records

if len(sys.argv) != 5:
	sys.exit("Usage: ion_plotGlobalCoverage.py <coverage.file> <bed.file> <output.folder> <analysis.name>");



coverageName = sys.argv[1];
targets = read_bed_records(sys.argv[2]);
outputDir = sys.argv[3];
analysisName = sys.argv[4];

x = [];
y = defaultdict(list);
for line in fileinput.input(coverageName):
	line = line.strip();
#	if not fileinput.isfirstline():
	tokens = line.split("\t");
	chr = tokens[0];
	loc = float(tokens[1]);
	depth = float(tokens[2]);
	depthF = float(tokens[3]);
	depthR = float(tokens[4]);
	for target in targets:
		if ((chr == targets[target].chrom) and (loc>=targets[target].chrStart) & (loc <= targets[target].chrStop)):
			targets[target].depth.append(depth);
			targets[target].depthF.append(depthF);
			targets[target].depthR.append(depthR);

fig = plt.figure(figsize=(15,15));
ax = fig.add_subplot(111);
offset = 0;
ticks = [];
labels = [];
for t in targets:
	y = targets[t].depth;
	if len(y)>0:
		ax.plot(arange(len(y)) + offset, y);
		ticks.append(offset + float(len(y))/2.0);
		offset = offset + len(y);
		labels.append(t);
ax.set_xticks(ticks);
ax.set_xticklabels(labels, size=10,rotation=90);
ax.grid();
plt.title(analysisName);
plt.savefig(outputDir + "/" + analysisName +"_global.png", format="PNG");


numTargets = float(len(targets));
numPlots = float(math.ceil(numTargets/60));
plotNum = 1;
i = 1;
numCols = 3;
numRows = math.ceil(numTargets/(numPlots*float(numCols)));

fig = plt.figure(figsize=(30,30));
plt.subplots_adjust(hspace=0.6);
for t in targets:
	if i > (numRows*numCols):
		plt.savefig(outputDir + "/" + analysisName + "_individual_" +str(plotNum)+".png", format="PNG");
		plotNum +=1;
		i = 1;
		plt.clf();
		plt.subplots_adjust(hspace=0.6);
	y = targets[t].depth;
	yF = targets[t].depthF;
	yR = targets[t].depthR;

	ax = fig.add_subplot(numRows,numCols,i);
	i+=1;
	#ax.bar(arange(len(y)), y, color='white');
	ax.bar(arange(len(y)), yF, color='red');
	ax.bar(arange(len(y)), yR, bottom=yF, color='green');
	ax.axis('tight');
	ax.grid();
	plt.title(t);
	#plt.ylabel("Start Position");
plt.savefig(outputDir + "/" + analysisName + "_individual_" +str(plotNum)+".png", format="PNG");



