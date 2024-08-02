#!/bin/bash
# Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved

write_file_links()
{
  local OUTDIR="$TSP_FILEPATH_PLUGIN_DIR"
  if [ -n ${1} ]; then
    OUTDIR=${1}
  fi
  local FILENAME="${HTML_RESULTS}"
  if [ -n ${2} ]; then
    FILENAME=${2}
  fi
  local OUTURL="$TSP_URLPATH_PLUGIN_DIR"
  if [ -n ${3} ]; then
    OUTURL=${3}
  fi
  local OUTFILE="${OUTDIR}/${FILENAME}"
  local HTSV="This is a tab-separated-values text file with a .xls filename extension."
  echo "<br/><div class='center' style='width:420px'>" >> "$OUTFILE"
  echo "<div class='grid-header' style='background-color:#DEF'><span style='float:left' class='ui-icon ui-icon-triangle-1-n'></span>" >> "$OUTFILE"
  echo "<span class='table-title'>File Links</span><span class='message'></span></div>" >> "$OUTFILE"

  # Quick links table
  echo "<div class='linkstable'>" >> "$OUTFILE"

  # Create IGV link
  echo "<div class='linkstable-row' id='IGV_LINK_ROW' row='even'>" >> "$OUTFILE"
  echo "<script type='text/javascript'>" >> "$OUTFILE"
  echo "  var locpath = window.location.pathname.substring(0,window.location.pathname.lastIndexOf('/'));" >> "$OUTFILE"
  echo "  var igvURL = window.location.protocol+'//'+window.location.host+'/auth'+locpath+'/igv.php3';" >> "$OUTFILE"
  echo "  var launchURL = window.location.protocol+'//'+window.location.host+'/IgvServlet/igv';" >> "$OUTFILE"
  echo "  var a = document.createElement('a');" >> "$OUTFILE"
  echo "  var linkText = document.createTextNode('View alignments in IGV.');" >> "$OUTFILE"
  echo "  a.appendChild(linkText);" >> "$OUTFILE"
  echo "  a.className = a.className+' flyhelp';" >> "$OUTFILE"
  echo "  a.title = 'Click to open the alignments and targets annotation in IGV.';" >> "$OUTFILE"
  echo "  a.href = launchURL + '?maxheap=1200M&locus=chr1&sessionURL='+igvURL;" >> "$OUTFILE"
  echo "  document.getElementById('IGV_LINK_ROW').appendChild(a);" >> "$OUTFILE"
  echo "</script>" >> "$OUTFILE"
  echo "</div>" >> "$OUTFILE"

  if [ -f "${OUTDIR}/$PLUGIN_OUT_STATSFILE" ]; then
    echo "<div class='linkstable-row' row='odd'>" >> "$OUTFILE"
    echo "<a class='flyhelp' href='${PLUGIN_OUT_STATSFILE}' title='Click to download a text file containing the summary statistics reported in the table above.'>Statistics summary file.</a><br/>" >> "$OUTFILE"
    echo "</div>" >> "$OUTFILE"
  fi
  if [ -f "${OUTDIR}/$PLUGIN_OUT_AMPLICAT" ]; then
    echo "<div class='linkstable-row' row='even'>" >> "$OUTFILE"
    echo "<a class='flyhelp' href='${PLUGIN_OUT_AMPLICAT}' title='Click to download a file containing read catagorization statistics for all individual amplicon targets.$HTSV'>Amplicons summary file.</a><br/>" >> "$OUTFILE"
    echo "</div>" >> "$OUTFILE"
  fi

  echo "</div></div><br/>"  >> "$OUTFILE"
}

