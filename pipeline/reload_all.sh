#!/bin/bash

./antismash_to_database.py
./clustercad_corrections.py
./clustercad_clean.py

mkdir -p /clusterCAD/pipeline/data/blast
cd /clusterCAD/pipeline
rm /clusterCAD/pipeline/data/blast/*
./blast_database.py
cd /clusterCAD/pipeline/data/blast
makeblastdb -in clustercad_subunits -parse_seqids -dbtype prot

cd /clusterCAD/pipeline
./process_aa_analysis.py
./generate_aa_plots.py
