#!/usr/bin/env python
import sys;
import fileinput;
from matplotlib import use
use("Agg");
import matplotlib.pyplot as plt;
import numpy as np;
import math;
import os;

if len(sys.argv) < 3:
    print("usage: ion_plotRepresentation.py <ampilcon.summary.file> <offtarget.summary.file> <output.dir>");
    exit(0);

repFile = sys.argv[1];
offFile = sys.argv[2];
outputDir = sys.argv[3];

labels = [];
ampgc = [];
offgc = [];

try:
	for line in fileinput.input(repFile):
		if not fileinput.isfirstline():
			tokens = line.strip().split(",");
                        labels.append(tokens[0]);
                        ampgc.append(float(tokens[15])*100.0);

except:
    print 'error while reading %s ' + repFile;
    exit(0);

try:
	for line in fileinput.input(offFile):
		if not fileinput.isfirstline():
			tokens = line.strip().split(",");
                        offgc.append(float(tokens[15])*100.0);

except:
    print 'error while reading %s ' + offFile;
    exit(0);


analysisName = os.path.basename(repFile);
analysisName = analysisName.split("_ampliconSummary.csv")[0];

gc = [ampgc, offgc];

fig = plt.figure(figsize=(20,10));
plt.boxplot(gc);
#plt.xticks(range(len(xticks)), xticks);
plt.xticks([1, 2], ['Amplicons', 'OffTargets']);
plt.grid();
plt.ylabel('GC%');
plt.xlabel('On/Off Targets');
plt.title('GC% of amplicons vs off-target regions');
plt.savefig(outputDir + "/" + analysisName + "_gc.on.vs.off.target.png", format="PNG");
