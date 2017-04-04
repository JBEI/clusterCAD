#!/bin/bash

function options {
cat << EOF
$0 options

This script processes a blastp query against the subuntis in the
ClusterCAD database.

Options:
 -q Query sequence
 -h Show help message
EOF
}

QUERY=
while getopts "q:h" OPTION; do
  case $OPTION in
    q)
      QUERY=$OPTARG;;
    h) 
      options
      exit 1;;
    ?) 
      options
      exit 1;;
  esac
done

blastp -db clustercad_subunits -max_target_seqs 10 -max_hsps 4 -outfmt "10 sseqid qstart qend sstart send evalue" -query <(echo -e ${QUERY})
