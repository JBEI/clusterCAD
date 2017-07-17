#!/bin/bash

./aggregate_final_structures.py
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

# loop generating plots until an error occurs, then stop
# script will exit with error code 1 if there are no
# plots left to generate
# this is a temporary hack to overcome a memory leak in the plot
# code
while [ $? -eq 0 ]; do
    ./generate_aa_plots.py
done
