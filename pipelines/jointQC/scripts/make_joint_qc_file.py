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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse input arguments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = "A simple script to demonstrate input arguments.")
parser.add_argument('--cellbender', type = str, help="")
parser.add_argument('--rna_metrics', type = str, help="")
parser.add_argument('--atac_metrics', type = str, help="")
parser.add_argument('--rna_whitelist', type = str, help="")
parser.add_argument('--atac_whitelist', type = str, help="")
parser.add_argument('--starsolo_dir', type = str, help="")
parser.add_argument('--outdir', type = str, help="")

args = parser.parse_args()

CELLBENDER = args.cellbender
RNA_METRICS = args.rna_metrics
ATAC_METRICS = args.atac_metrics
RNA_BARCODE_WHITELIST = args.rna_whitelist
ATAC_BARCODE_WHITELIST = args.atac_whitelist
GENE_FULL_EXON_OVER_INTRON_COUNTS = args.starsolo_dir + 'GeneFull_ExonOverIntron/raw'
GENE_COUNTS = args.starsolo_dir + 'Gene/raw'
OUTDIR = args.outdir

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in QC output from both modalities
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

qc_df = metrics

qc_df.to_csv(OUTDIR + "joint_qc.csv", index=False)