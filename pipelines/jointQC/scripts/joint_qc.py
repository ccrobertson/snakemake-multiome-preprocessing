from cellbender.remove_background.downstream import load_anndata_from_input_and_output
from cellbender.remove_background.downstream import anndata_from_h5
import scanpy as sc
import pandas as pd
from scipy.io import mmread
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import upsetplot

# import tables
# import anndata
# from typing import Dict, Optional
# 
# import scipy.sparse as sp
# from scipy import io
# import glob



### DEFINE INPUT FILES
ROOT = '/scratch/scjp_root/scjp0/ccrober/preprocessing/9240-VD'
CELLBENDER = ROOT+'/'+'work/multiome-rna/results/cellbender/9240-VD-1-hg38.cellbender_FPR_0.1.h5'
RNA_METRICS = ROOT+'/'+'work/multiome-rna/results/qc/9240-VD-1-hg38.qc.txt'
ATAC_METRICS = ROOT+'/'+'work/multiome-atac/results/ataqv/single-nucleus/9240-VD-1-hg38.txt'
RNA_BARCODE_WHITELIST = ROOT+'/'+'pipelines/snRNAseq-NextFlow/737K-arc-v1.txt.gz'
ATAC_BARCODE_WHITELIST = ROOT+'/'+'pipelines/snATACseq-NextFlow/737K-arc-v1.txt'
GENE_FULL_EXON_OVER_INTRON_COUNTS = ROOT+'/'+'work/multiome-rna/results/starsolo/9240-VD-1-hg38/9240-VD-1-hg38.Solo.out/GeneFull_ExonOverIntron/raw'
GENE_COUNTS = ROOT+'/'+'work/multiome-rna/results/starsolo/9240-VD-1-hg38/9240-VD-1-hg38.Solo.out/Gene/raw'

### DEFINE THRESHOLDS
THRESHOLD_RNA_MIN_UMI = 150
THRESHOLD_RNA_MAX_MITO = 0.1
THRESHOLD_ATAC_MIN_HQAA = 5000
THRESHOLD_ATAC_MIN_TSS_ENRICHMENT = 2.5
THRESHOLD_CELLBENDER_MIN_CELL_PROBABILITY = 0.99
THRESHOLD_POST_CELLBENDER_UMIS = 100
THRESHOLD_EXON_GENE_BODY_RATIO = 0.5


### FUNCTIONS TO EXTRA INFO            
def cellbender_anndata_to_cell_probability(a):
    return a.obs.cell_probability

def cellbender_anndata_to_sparse_matrix(adata, min_cell_probability=0):
    barcodes = adata.obs[adata.obs.cell_probability>=min_cell_probability].index.to_list()
    features = adata.var.gene_id.to_list()
    matrix = adata[adata.obs.cell_probability>=min_cell_probability].X.transpose()
    return {'features': features, 'barcodes': barcodes, 'matrix': matrix}

def umi_count_after_decontamination(adata):
    x = cellbender_anndata_to_sparse_matrix(adata)
    return dict(zip(x['barcodes'], x['matrix'].sum(axis=0).tolist()[0]))


### PLOTTING FUNCTIONS
def barcode_rank_plot(metrics, ax):
    df = metrics.sort_values('rna_umis', ascending=False)
    df['barcode_rank'] = range(1, len(df) + 1)
    sns.scatterplot(x='barcode_rank', y='rna_umis', data=df, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Barcode rank')
    ax.set_ylabel('UMIs')
    return ax

def rna_umis_vs_rna_mito_plot(metrics, ax):
    sns.scatterplot(x='rna_umis', y='rna_fraction_mitochondrial', data=metrics, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.02, s=3)
    ax.set_xscale('log')
    ax.set_xlabel('UMIs')
    ax.set_ylabel('Fraction mito. (RNA)')
    return ax

def rna_umis_vs_exon_to_full_gene_body_ratio(metrics, ax):
    sns.scatterplot(x='rna_umis', y='rna_exon_to_full_gene_body_ratio', data=metrics, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.02, s=3)
    ax.set_xscale('log')
    ax.set_xlabel('UMIs')
    ax.set_ylabel('Exon/full-gene-body ratio (RNA)')
    return ax

def cellbender_fraction_removed(metrics, ax):
    sns.scatterplot(x='rna_umis', y='fraction_cellbender_removed', data=metrics, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.05)
    ax.set_xscale('log')
    ax.set_xlabel('UMIs')
    ax.set_ylabel('Fraction ambient')
    return ax

def cellbender_cell_probabilities(metrics, ax):
    sns.histplot(x='cell_probability', data=metrics[(metrics.filter_rna_min_umi) & (metrics.filter_rna_max_mito)], ax=ax, bins=20)
    ax.set_xlabel('Cellbender cell prob.\nfor nuclei passing UMI and mito. thresholds')
    return ax

def rna_umis_vs_atac_hqaa_plot(metrics, ax):
    sns.scatterplot(x='rna_umis', y='atac_hqaa', data=metrics, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.02, s=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('UMIs (RNA)')
    ax.set_ylabel('Pass filter reads (ATAC)')
    return ax

def atac_hqaa_vs_atac_tss_enrichment_plot(metrics, ax):
    sns.scatterplot(x='atac_hqaa', y='atac_tss_enrichment', data=metrics, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.02, s=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Pass filter reads (ATAC)')
    ax.set_ylabel('TSS enrichment')
    return ax

### BUILD BARCODE MAP
# ATAC --> RNA barcode mappings
rna_barcodes = pd.read_csv(RNA_BARCODE_WHITELIST, header=None)[0].to_list()
atac_barcodes = pd.read_csv(ATAC_BARCODE_WHITELIST, header=None)[0].to_list()
atac_to_rna = dict(zip(atac_barcodes, rna_barcodes))

### Calculate ratio of exonic vs full gene body reads
# exons only
gene_mat = mmread(os.path.join(GENE_COUNTS, 'matrix.mtx'))
gene_umis_per_barcode = gene_mat.sum(axis=0).tolist()[0]

# includes introns
gene_full_mat = mmread(os.path.join(GENE_FULL_EXON_OVER_INTRON_COUNTS, 'matrix.mtx'))
gene_full_umis_per_barcode = gene_full_mat.sum(axis=0).tolist()[0]

barcodes = pd.read_csv(os.path.join(GENE_COUNTS, 'barcodes.tsv'), header=None)[0]
assert(all(barcodes == pd.read_csv(os.path.join(GENE_FULL_EXON_OVER_INTRON_COUNTS, 'barcodes.tsv'), header=None)[0]))

exon_to_full_gene_body_ratio = pd.DataFrame({'barcode': barcodes, 'gene': gene_umis_per_barcode, 'gene_full': gene_full_umis_per_barcode})
exon_to_full_gene_body_ratio['exon_to_full_gene_body_ratio'] = exon_to_full_gene_body_ratio.gene / exon_to_full_gene_body_ratio.gene_full
exon_to_full_gene_body_ratio.head()

adata = anndata_from_h5(CELLBENDER, analyzed_barcodes_only=True)
rna_metrics = pd.read_csv(RNA_METRICS, sep='\t')
rna_metrics = rna_metrics[rna_metrics.barcode!='-']
rna_metrics = rna_metrics.merge(exon_to_full_gene_body_ratio)

atac_metrics = pd.read_csv(ATAC_METRICS, sep='\t', header=None, names=['barcode', 'metric', 'value'])
KEEP_ATAC_METRICS = ['median_fragment_length', 'hqaa', 'max_fraction_reads_from_single_autosome', 'percent_mitochondrial', 'tss_enrichment']
atac_metrics = atac_metrics[atac_metrics.metric.isin(KEEP_ATAC_METRICS)].pivot(index='barcode', columns='metric', values='value')
atac_metrics.hqaa = atac_metrics.hqaa.astype(int)
atac_metrics.max_fraction_reads_from_single_autosome = atac_metrics.max_fraction_reads_from_single_autosome.astype(float)
atac_metrics.median_fragment_length = atac_metrics.median_fragment_length.astype(float)
atac_metrics.percent_mitochondrial = atac_metrics.percent_mitochondrial.astype(float)
atac_metrics.tss_enrichment = atac_metrics.tss_enrichment.astype(float)
atac_metrics['fraction_mitochondrial'] = atac_metrics.percent_mitochondrial / 100

atac_metrics.index = atac_metrics.index.map(atac_to_rna)
metrics = rna_metrics.set_index('barcode').rename(columns=lambda x: 'rna_' + x).join(atac_metrics.rename(columns=lambda x: 'atac_' + x))

metrics = metrics.reset_index()
cell_probability = cellbender_anndata_to_cell_probability(adata)
post_cellbender_umis = umi_count_after_decontamination(adata)

metrics['cell_probability'] = metrics.barcode.map(lambda x: cell_probability[x] if x in cell_probability else np.nan)
metrics['post_cellbender_umis'] = metrics.barcode.map(lambda x: post_cellbender_umis[x] if x in post_cellbender_umis else np.nan)
metrics['fraction_cellbender_removed'] = (metrics.rna_umis - metrics.post_cellbender_umis) / metrics.rna_umis
metrics.head()

# apply QC thresholds
metrics['filter_cellbender_cell_probability'] = metrics.cell_probability >= THRESHOLD_CELLBENDER_MIN_CELL_PROBABILITY
metrics['filter_rna_min_umi'] = metrics.rna_umis >= THRESHOLD_RNA_MIN_UMI
metrics['filter_rna_postcellbender_min_umi'] = metrics.post_cellbender_umis >= THRESHOLD_POST_CELLBENDER_UMIS
metrics['filter_rna_exon_to_full_gene_body_ratio'] = metrics.rna_exon_to_full_gene_body_ratio <= THRESHOLD_EXON_GENE_BODY_RATIO
metrics['filter_rna_max_mito'] = metrics.rna_fraction_mitochondrial <= THRESHOLD_RNA_MAX_MITO
metrics['filter_atac_min_hqaa'] = metrics.atac_hqaa >= THRESHOLD_ATAC_MIN_HQAA
metrics['filter_atac_min_tss_enrichment'] = metrics.atac_tss_enrichment >= THRESHOLD_ATAC_MIN_TSS_ENRICHMENT

metrics['pass_all_filters'] = metrics.filter(like='filter_').all(axis=1)

# List of pass-QC barcodes
pass_qc_nuclei = list(sorted(metrics[metrics.pass_all_filters].barcode.to_list()))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot QC metrics 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig, axs = plt.subplots(ncols=7, figsize=(7*5, 4))

ax = axs[0]
barcode_rank_plot(metrics, ax)
ax.axhline(THRESHOLD_RNA_MIN_UMI, color='red', ls='--')

ax = axs[1]
rna_umis_vs_rna_mito_plot(metrics, ax)
ax.axhline(THRESHOLD_RNA_MAX_MITO, color='red', ls='--')
ax.axvline(THRESHOLD_RNA_MIN_UMI, color='red', ls='--')

ax = axs[2]
rna_umis_vs_exon_to_full_gene_body_ratio(metrics, ax)
ax.axhline(THRESHOLD_EXON_GENE_BODY_RATIO, color='red', ls='--')
ax.axvline(THRESHOLD_RNA_MIN_UMI, color='red', ls='--')
ax.set_xlim(left=0.8*THRESHOLD_RNA_MIN_UMI)

ax = axs[3]
cellbender_fraction_removed(metrics, ax)

ax = axs[4]
cellbender_cell_probabilities(metrics, ax)

ax = axs[5]
rna_umis_vs_atac_hqaa_plot(metrics, ax)
ax.axhline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--')
ax.axvline(THRESHOLD_RNA_MIN_UMI, color='red', ls='--')

ax = axs[6]
atac_hqaa_vs_atac_tss_enrichment_plot(metrics, ax)
ax.axvline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--')
ax.axhline(THRESHOLD_ATAC_MIN_TSS_ENRICHMENT, color='red', ls='--')

fig.suptitle('{:,} pass QC nuclei'.format(len(pass_qc_nuclei)))
fig.tight_layout()
plt.savefig('work/joint_qc/joint_qc_plots.png', format='jpeg', dpi=300)
    
    # Plot the number of nuclei passing each filter
    #fig, ax = plt.subplots(figsize=(7, 6))
    #ax.remove()

    
    # for_upset = metrics.filter(like='filter_').rename(columns=lambda x: 'pass_' + x)
    # for_upset = for_upset.groupby(for_upset.columns.to_list()).size()
    # upsetplot.plot(for_upset, fig=fig, sort_by='cardinality', show_counts=True)


with PdfPages('work/joint_qc/output_plots_test.pdf') as pdf:

    fig, axs = plt.subplots(ncols=3)
    axs[0].plot([1, 2, 3, 4], [1, 4, 2, 3])
    axs[1].plot([1, 2, 3, 4], [4, 5, 6, 7])
    pdf.savefig()



fig, ax = plt.subplots()
barcode_rank_plot(metrics, ax)
ax.axhline(THRESHOLD_RNA_MIN_UMI, color='red', ls='--')
plt.savefig('work/joint_qc/barcode_rank_plot.png', format='jpeg', dpi=300)

