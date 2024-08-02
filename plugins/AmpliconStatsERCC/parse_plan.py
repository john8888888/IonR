#!/usr/bin/env python
# Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved
import os
import sys
import json

   
if __name__ == '__main__':
  try:
    JSON_INPUT = json.load( open(sys.argv[1], "r") )
    PLAN_INFO = JSON_INPUT['plan']
    PLUGIN_CONFIG = JSON_INPUT['pluginconfig']
  except:
    print ";;;"
    sys.exit(0) 

  if PLAN_INFO:
     runtype = PLAN_INFO['runType']
     regionf = PLAN_INFO['bedfile']

     analysistype = 'Restart'
     if 'analysis_type' in PLUGIN_CONFIG:
       analysistype = PLUGIN_CONFIG['analysis_type']

     offtarget = 'Off'
     if 'offtarget_analysis' in PLUGIN_CONFIG:
       offtarget = PLUGIN_CONFIG['offtarget_analysis']

     prenumber = ''
     if 'number' in PLUGIN_CONFIG:
       prenumber = PLUGIN_CONFIG['number']
     if ( runtype == 'AMPS' ):
        runtype = "ampliseq"       
     elif ( runtype == 'TARS' ):
        runtype = "targetseq"
     elif ( runtype == 'WGNM' ):
        runtype = "fullgenome"
     else:
        runtype = "wrongtype"  

     #print "%s;%s;%s;%s" % (runtype, regionf, analysistype, offtarget)
     print "%s;%s;%s;%s;%s" % (runtype, regionf, analysistype, offtarget, prenumber)
  else:
     print ";;;;"
