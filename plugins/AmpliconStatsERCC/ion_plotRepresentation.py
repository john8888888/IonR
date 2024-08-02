#!/usr/bin/env python
import sys;
import fileinput;
from matplotlib import use
use("Agg");
import matplotlib.pyplot as plt;
import numpy as np;
import math;
import os;
import datetime;
from subprocess import call;

#import Gnuplot, Gnuplot.funcutils;
#from subprocess import call;

if len(sys.argv) < 2:
    print("usage: ion_plotRepresentation.py <ampilcon.summary.file> <output.dir> <num.amplicons.per.chart>");
    exit(0);

repFile = sys.argv[1];
outputDir = sys.argv[2];
plugindir = sys.argv[3];

#if len(sys.argv) == 4:
#	maxTargets = int(sys.argv[3]);
#else:
#maxTargets = 100;



fwd = [];
rev = [];
labels = [];

xd = [];
yd = [];

now = datetime.datetime.now()
print "Beginning to process amplicon summary: "
print now.strftime("%Y-%m-%d %H:%M:%S")

try:
	for line in fileinput.input(repFile):
		if not fileinput.isfirstline():
			tokens = line.strip().split(",");
			labels.append(tokens[0]);
			fwd.append(float(tokens[1]));
			rev.append(float(tokens[2]));
                        xd.append(float(tokens[15])*100.0);
                        yd.append(float(tokens[11]));

except:
    print 'error while reading %s ' + repFile;
    exit(0);

now = datetime.datetime.now()
print "finished processing amplicon summary: "
print now.strftime("%Y-%m-%d %H:%M:%S")

analysisName = os.path.basename(repFile);
analysisName = analysisName.split("_ampliconSummary.csv")[0];
print "analysis name is " + analysisName;

repPNG = outputDir + "/" + analysisName + "_representation.png"
command = "Rscript " + plugindir + "/ion_plotRepresentation.R " + repFile + " " + repPNG;
print command;
call(command, shell=True);

numTargets = float(len(fwd));
#numPlots = math.ceil(numTargets/maxTargets);
#numTargetsPerPlot = math.ceil(numTargets/numPlots);

fwd = (np.array(fwd));
rev = (np.array(rev));
total = fwd + rev;
ind = np.argsort(total);

labels = np.array(labels);
total = total[ind];
fwd = fwd[ind];
rev = rev[ind];
labels = labels[ind];

now = datetime.datetime.now()
print "Beginning to get ci95: "
print now.strftime("%Y-%m-%d %H:%M:%S")

#maxVal = np.max(fwd+rev);
mu = np.mean(total);

maxVal = max(np.max(fwd), np.max(rev));
muFwd = np.mean(fwd);
muRev = np.mean(rev);

width = 1;
negRev = 0 - rev;
x = np.arange(len(fwd));

now = datetime.datetime.now()
print "Beginning to plot: "
print now.strftime("%Y-%m-%d %H:%M:%S")


fig = plt.figure(figsize=(20,10));
maxLog = max(np.max(np.log10(fwd)), np.max(np.log10(rev)));
plt.bar(x, np.log10(fwd), width, edgecolor='orange', color='orange', alpha=0.7);
plt.bar(x, 0-np.log10(rev), width, edgecolor='green', color='green', alpha=0.7);

#plt.yscale('log');
plt.title('Total Number of Amplicon Reads (sorted)\n Log Scale');
plt.ylim([-maxLog, maxLog]);
plt.xlim([0, len(fwd)]);
plt.savefig(outputDir + "/" + analysisName + "_log_representation.png", format="PNG");

now = datetime.datetime.now()
print "mat pyplot finishes log representation plot: "
print now.strftime("%Y-%m-%d %H:%M:%S")


den = (fwd + rev);
repBias = fwd / den;
repBias[np.isinf(repBias)] = 0;
repBias[np.isnan(repBias)] = 0;
repBias = np.sort(repBias);
numPassRepBias = float(sum((repBias>=0.3) & (repBias<=0.7)));
percPassRepBias = numPassRepBias/float(len(repBias));

plt.clf();
plt.axhline(y=0.3, color='orange', alpha=0.2, linewidth=4)
plt.axhline(y=0.7, color='orange', alpha=0.2, linewidth=4)
plt.axhline(y=0.4, color='g', alpha=0.2, linewidth=4)
plt.axhline(y=0.6, color='g', alpha=0.2, linewidth=4)
plt.plot(repBias, 'o', color='red');
plt.grid();
plt.ylim([0, 1]);
plt.xlim([0, len(fwd)]);
plt.ylabel(r'Probability of forward read $(\beta)$');
plt.xlabel('Sorted amplicon');
plt.title(('Plot Showing Representation Bias by Amplicon\n%2.2f%% amplicons pass [%s], %s')%(100.*percPassRepBias, r'$0.3 \leq \beta \leq 0.7$', r'$\beta=\frac{FWD}{FWD+REV}$'));
plt.savefig(outputDir + "/" + analysisName + "_repBias.png", format="PNG");

now = datetime.datetime.now()
print "mat pyplot finished bias plot: "
print now.strftime("%Y-%m-%d %H:%M:%S")


#plot sucFailvsGC.png
outFile = analysisName + '_sucFailvsGC.png';

fig = plt.figure(figsize=(20,10));
title =  'Number of amplicon success(green)/Fail(red) to pass 0.2x mean reads per base vs GC content';
plt.title(title);
plt.grid();
plt.xlabel('GC percentage');
plt.ylabel('Number of success/fail');
newd = dict();
llimit = int(round(min(xd)));
hlimit = int(round(max(xd)));
bins = np.linspace(llimit, hlimit, hlimit - llimit + 1);
for i in range(len(xd)):
    d1 = str(int(round(xd[i])));
    newd[d1] = [0, 0];
for i in range(len(xd)):
    d1 = str(int(round(xd[i])));
    d2 = str(int(yd[i]));
    if d2 == '0':
        newd[d1][0] += 1;
        newd[d1][1] += 0;
    elif d2 == '1':
        newd[d1][1] += 1;
        newd[d1][0] += 0;
    else:
        print "wrong value of " + d2;
        exit(0);

newx = [];
newy0 = [];
newy1 = [];
for k in sorted(newd):
    newx.append(int(k));
    newy0.append(newd[k][0]);
    newy1.append(newd[k][1]);

newy0Normalized = [];
newy1Normalized = [];
for k in sorted(newd):
    norm = float(newd[k][0] +  newd[k][1])/100.;
    newy0Normalized.append(float(newd[k][0])/norm);
    newy1Normalized.append(float(newd[k][1])/norm);

now = datetime.datetime.now()
print "mat pyplot begins to draw sucFailvsGC plot: "
print now.strftime("%Y-%m-%d %H:%M:%S")

#presentation format 1
plt.bar(newx, newy0, bottom=0, color='red');
plt.bar(newx, newy1, bottom=newy0, color='green');
plt.savefig(outFile, format="PNG");

#plot sucFailvsGC.png
outFile = analysisName + '_sucFailRatevsGC.png';

fig = plt.figure(figsize=(20,10));
title =  'Percentage of amplicon success(green)/Fail(red) to pass 0.2x mean reads per base vs GC content';
plt.title(title);
plt.grid();
plt.xlabel('GC percentage');
plt.ylabel('Percentage of success/fail');
#presentation format 3
plt.plot(newx, newy0Normalized, color='red', linewidth=2, marker='1');
plt.plot(newx, newy1Normalized, color='green', linewidth=2, marker='2');
plt.savefig(outFile, format="PNG");

now = datetime.datetime.now()
print "mat pyplot finishes sucFailvsGC plot: "
print now.strftime("%Y-%m-%d %H:%M:%S")

