#!/bin/bash

function options {
cat << EOF
$0 options

This script downloads and processes the MIBiG files to populate the
ClusterCAD database.

Options:
 -v Version of MiBiG to download
 -h Show help message
EOF
}

VERSION=
while getopts "v:h" OPTION; do
  case $OPTION in
    v)
      VERSION=$OPTARG;;
    h) 
      options
      exit 1;;
    ?) 
      options
      exit 1;;
  esac
done

#######################
# Download MIBiG data #
#######################

mkdir -p ./data/mibig/raw

wget -O ./data/mibig/mibig_json_${VERSION}.tar.gz http://dl.secondarymetabolites.org/mibig/mibig_json_${VERSION}.tar.gz
wget -O ./data/mibig/mibig_gbk_${VERSION}.tar.gz http://dl.secondarymetabolites.org/mibig/mibig_gbk_${VERSION}.tar.gz

tar -xvf ./data/mibig/mibig_json_${VERSION}.tar.gz -C ./data/mibig/raw
tar -xvf ./data/mibig/mibig_gbk_${VERSION}.tar.gz -C ./data/mibig/raw

#################################
# Process antiSMASH annotations #
#################################

mkdir -p ./data/antismash/raw
# Generate antismash output in folder ./data/antismash/raw (previously mibig.tar.gz)

mkdir -p ./data/antismash/split

#awk 'BEGIN{ RS="//\n"; ORS="//"; } {split($2, array, ";"); fname=array[1]; print > "./data/antismash/split/" fname ".embl" }' ./data/antismash/raw/BGC0000001.1.final.embl

