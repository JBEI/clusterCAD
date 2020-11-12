import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import os,sys
import glob
antismash_path = "antismash_db"
filelist = glob.glob(os.path.join(antismash_path, '*.gbk'))
#print(filelist)
total_count = []
for file in filelist:
    record = SeqIO.read(file, "genbank")
    count = 0
    for feature in record.features:
        if not feature.type == 'CDS' or 'sec_met' not in feature.qualifiers.keys():
            continue
        if 't1pks' not in feature.qualifiers['sec_met'][0]:
            continue
        count+= 1
    print('added count', count, 'for', record.id)
    total_count.append(count)

print('average subunits per pks is', np.mean(total_count), 'and the median is', np.median(total_count))
