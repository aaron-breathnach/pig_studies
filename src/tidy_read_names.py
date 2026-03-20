import glob
import os

filenames = glob.glob('reads/raw/*_1.fastq.gz')

sample_ids = [os.path.basename(x.replace('_1.fastq.gz', '')) for x in glob.glob('reads/raw/*_1.fastq.gz')]

for sample_id in sample_ids:
    for r in ['1', '2']:
       old = f'reads/raw/{sample_id}_{r}.fastq.gz'
       new = f'reads/raw/{sample_id}_R{r}_001.fastq.gz'
       cmd = f'mv {old} {new}'
       os.system(cmd)
