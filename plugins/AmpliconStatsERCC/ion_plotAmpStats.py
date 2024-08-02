#!/usr/bin/env python
import sys;
import fileinput;
from matplotlib import use
use("Agg");
import matplotlib.pyplot as plt;
import numpy as np;
import math;
import os;

if len(sys.argv) != 3:
    print("usage: ion_plotAmpStats.py <ampilcon.summary.file> <output.dir>");
    exit(0);

repFile = sys.argv[1];
outputDir = sys.argv[2];

fwd = [];
rev = [];
gc = [];
length = [];
labels = [];
try:
	for line in fileinput.input(repFile):
		if not fileinput.isfirstline():
			tokens = line.strip().split(",");
			name = tokens[0];
			labels.append(name);
			fwd.append(float(tokens[1]));
			rev.append(float(tokens[2]));
			gc.append(float(tokens[15]));
			length.append(float(len(tokens[16])));
except:
	exit(0);

analysisName = os.path.basename(repFile);
analysisName = analysisName.split("_ampliconSummary.csv")[0];


fig = plt.figure(figsize=(20,10));
plt.plot(fwd, gc, 'o', color='lightblue', alpha=0.5);
plt.ylim([0, 1]);
plt.grid();
plt.title('Forward Reads\nMean Amplicon GC = %2.2f, Max. Amplicon GC = %2.2f, Min. Amplicon GC = %2.2f' % (np.mean(gc), np.max(gc), np.min(gc)));
plt.xlabel('Number of Forward Reads');
plt.ylabel('Amplicon GC');
plt.savefig(outputDir + "/" + analysisName + "_gc_fwd.png", format="PNG");

plt.clf();
plt.plot(rev, gc, 'o', color='g', alpha=0.5);
plt.ylim([0, 1]);
plt.grid();
plt.title('Reverse Reads\nMean Amplicon GC = %2.2f, Max. Amplicon GC = %2.2f, Min. Amplicon GC = %2.2f' % (np.mean(gc), np.max(gc), np.min(gc)));
plt.xlabel('Number of Reverse Reads');
plt.ylabel('Amplicon GC');
plt.savefig(outputDir + "/" + analysisName + "_gc_rev.png", format="PNG");

plt.clf();
plt.plot(fwd, length, 'o', color='orange', alpha=0.5);
plt.grid();
plt.title('Forward Reads\nMean Amplicon Length = %2.2f, Max. Amplicon Length = %2.2f, Min. Amplicon Length = %2.2f' % (np.mean(length), np.max(length), np.min(length)));
plt.xlabel('Number of Forward Reads');
plt.ylabel('Amplicon Length');
plt.savefig(outputDir + "/" + analysisName + "_len_fwd.png", format="PNG");

plt.clf();
plt.plot(rev, length, 'o', color='purple', alpha=0.5);
plt.grid();
plt.title('Reverse Reads\nMean Amplicon Length = %2.2f, Max. Amplicon Length = %2.2f, Min. Amplicon Length = %2.2f' % (np.mean(length), np.max(length), np.min(length)));
plt.xlabel('Number of Reverse Reads');
plt.ylabel('Amplicon Length');
plt.savefig(outputDir + "/" + analysisName + "_len_rev.png", format="PNG");
