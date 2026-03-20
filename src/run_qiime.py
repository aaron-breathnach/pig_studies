import argparse
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--fastq', help='The directory containing the paired reads', required=True)
parser.add_argument('--qiime', help='The output directory', required=True)
parser.add_argument('--metadata', help='Sample metadata', required=True)
parser.add_argument('--threads', help='Number of threads to use', default=1)
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

step_01 = 'Rscript src/make_manifest.R {fastq} {qiime}'

step_02 = '''
qiime tools import 
--type 'SampleData[PairedEndSequencesWithQuality]' 
--input-format PairedEndFastqManifestPhred33V2 
--input-path {qiime}/manifest.tsv 
--output-path {qiime}/demux_seqs.qza
'''.replace('\n', '')

step_03 = '''
qiime cutadapt trim-paired 
--i-demultiplexed-sequences {qiime}/demux_seqs.qza 
--p-cores {threads} 
--p-front-f CCTACGGGNGGCWGCAG 
--p-front-r GACTACHVGGGTATCTAATCC 
--o-trimmed-sequences {qiime}/demux_seqs.trimmed.qza
'''.replace('\n', '')

step_04 = 'qiime demux summarize --i-data {qiime}/demux_seqs.trimmed.qza --o-visualization {qiime}/demux_seqs.trimmed.qzv'

step_05 = 'qiime tools export --input-path {qiime}/demux_seqs.trimmed.qzv --output-path {qiime}/demultiplex_summary'

step_06 = 'Rscript src/get_cutoffs.R {qiime}/demultiplex_summary {qiime}'

step_07 = '''
qiime dada2 denoise-paired 
--i-demultiplexed-seqs {qiime}/demux_seqs.trimmed.qza 
--p-trim-left-f 0 
--p-trim-left-r 0 
--p-trunc-len-f {trunc_len_f} 
--p-trunc-len-r {trunc_len_r} 
--p-n-threads 0 
--o-table {qiime}/table.qza 
--o-representative-sequences {qiime}/rep-seqs.qza 
--o-denoising-stats {qiime}/denoising-stats.qza 
--o-base-transition-stats {qiime}/base-transition-stats.qza
'''.replace('\n', '')

step_08 = '''
qiime feature-table summarize 
--i-table {qiime}/table.qza 
--o-visualization {qiime}/table.qzv 
--m-sample-metadata-file {metadata}
'''.replace('\n', '')

step_09 = '''
qiime feature-table tabulate-seqs 
--i-data {qiime}/rep-seqs.qza 
--o-visualization {qiime}/rep-seqs.qzv
'''.replace('\n', '')

step_10 = '''
qiime metadata tabulate 
--m-input-file {qiime}/denoising-stats.qza 
--o-visualization {qiime}/denoising-stats.qzv
'''.replace('\n', '')

step_11 = '''
qiime phylogeny align-to-tree-mafft-fasttree 
--i-sequences {qiime}/rep-seqs.qza 
--o-alignment {qiime}/aligned-rep-seqs.qza 
--o-masked-alignment {qiime}/masked-aligned-rep-seqs.qza 
--o-tree {qiime}/unrooted-tree.qza 
--o-rooted-tree {qiime}/rooted-tree.qza
'''.replace('\n', '')

step_12 = '''
qiime feature-classifier classify-sklearn 
--i-classifier {qiime}/gg-13-8-99-515-806-nb-classifier.qza 
--i-reads {qiime}/rep-seqs.qza 
--o-classification {qiime}/taxonomy.qza
'''.replace('\n', '')

step_13 = 'qiime tools export --input-path {qiime}/denoising-stats.qza --output-path {qiime}'

step_14 = 'Rscript src/get_p_sampling_depth.R {qiime}/stats.tsv {qiime}'

step_15 = '''
qiime diversity core-metrics-phylogenetic 
--i-phylogeny {qiime}/rooted-tree.qza 
--i-table {qiime}/table.qza 
--p-sampling-depth {p_sampling_depth} 
--m-metadata-file {metadata} 
--output-dir {qiime}/core-metrics-results
'''.replace('\n', '')

def run_qiime(fastq, qiime, metadata, threads, rerun=False):
    dictionary = {
        'step01': {'cmd': step_01.format(fastq=fastq, qiime=qiime), 'target': '{}/manifest.tsv'.format(qiime)},
        'step02': {'cmd': step_02.format(qiime=qiime), 'target': '{}/demux_seqs.qza'.format(qiime)},
        'step03': {'cmd': step_03.format(qiime=qiime, threads=threads), 'target': '{}/demux_seqs.trimmed.qza'.format(qiime)},
        'step04': {'cmd': step_04.format(qiime=qiime), 'target': '{}/demux_seqs.trimmed.qzv'.format(qiime)},
        'step05': {'cmd': step_05.format(qiime=qiime), 'target': '{}/demultiplex_summary/index.html'.format(qiime)},
        'step06': {'cmd': step_06.format(qiime=qiime), 'target': '{}/trunc_len.txt'.format(qiime)},
        'step07': {'cmd': step_07, 'target': '{}/rep-seqs.qza'.format(qiime)},
        'step08': {'cmd': step_08.format(qiime=qiime, metadata=metadata), 'target': '{}/table.qzv'.format(qiime)},
        'step09': {'cmd': step_09.format(qiime=qiime), 'target': '{}/rep-seqs.qzv'.format(qiime)},
        'step10': {'cmd': step_10.format(qiime=qiime), 'target': '{}/denoising-stats.qzv'.format(qiime)},
        'step11': {'cmd': step_11.format(qiime=qiime), 'target': '{}/aligned-rep-seqs.qza'.format(qiime)},
        'step12': {'cmd': step_12.format(qiime=qiime), 'target': '{}/taxonomy.qza'.format(qiime)},
        'step13': {'cmd': step_13.format(qiime=qiime), 'target': '{}/stats.tsv'.format(qiime)},
        'step14': {'cmd': step_14.format(qiime=qiime), 'target': '{}/p_sampling_depth.txt'.format(qiime)},
        'step15': {'cmd': step_15, 'target': '{}/core-metrics-results/shannon_vector.qza'.format(qiime)}
    }
    cmds = []
    for i in dictionary:
        if not os.path.exists(dictionary[i]['target']) or rerun:
            cmds.append(dictionary[i]['cmd'])
        else:
            cmds.append(None)
    return(cmds)

if __name__ == '__main__':
    greengenes = '{}/gg-13-8-99-515-806-nb-classifier.qza'.format(args.qiime)
    if not os.path.exists(greengenes):
        url = 'https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza'
        wget = 'wget {} -O {}'.format(url, greengenes)
        subprocess.run(wget, shell=True)
    steps = run_qiime(args.fastq, args.qiime, args.metadata, args.threads, args.rerun)
    for step in steps[0:5]:
        if step != None:
            os.system(step)
    trunc_len = [int(x) for x in open('{}/trunc_len.txt'.format(args.qiime), 'r').readlines()]
    step_07 = steps[6].format(qiime=args.qiime, trunc_len_f=trunc_len[0], trunc_len_r=trunc_len[1])
    os.system(step_07)
    for step in steps[7:13]:
        os.system(step)
    p_sampling_depth = int(open('{}/p_sampling_depth.txt'.format(args.qiime), 'r').readlines()[0])
    step_15 = steps[14].format(qiime=args.qiime, p_sampling_depth=p_sampling_depth-1, metadata=args.metadata)
    os.system(step_15)
