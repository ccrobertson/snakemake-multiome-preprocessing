## Description
This pipeline performs 
1. Automated detection of nuclei-containing droplets --> barcode_keeplist
2. Joint qc plots to assess quality of droplets in barcode_keeplist

As input this requires a config table with the following columns:
- library
- path_type
- path

Where the following path types are provided:
- starsolo_dir
- cellbender_h5
- rna_metrics
- atac_metrics
- rna_whitelist
- atac_whitelist

## Run the pipeline
```
bash run.sh
```


## For conciseness, will define models as
* Model 1 - RNA umis + RNA mitochondrial fraction
* Model 2 - RNA umis + ATAC HQAA
* Model 3 - RNA umis + ATAC TSS enrichment
* Model 4 - RNA umis + RNA mitochondrial fraction + ATAC TSS enrichment
* Model 5 - RNA umis + RNA exonic fraction

## Compile the images for sharing
```
module load Bioinformatics  gcc/10.3.0-k2osx5y
module load imagemagick/7.1.0-33-xlaibek

model_1="gmm__rna_umis__rna_fraction_mitochondrial"
model_2="gmm__rna_umis__atac_tss_enrichment"
model_3="gmm__rna_umis__rna_exon_to_full_gene_body_ratio"
model_4="gmm__rna_umis__rna_fraction_mitochondrial__atac_tss_enrichment"
model_5="gmm__rna_umis__rna_fraction_mitochondrial__rna_exon_to_full_gene_body_ratio"

convert -limit memory 64 results/joint_qc_gmm/*/*${model_1}_bestfit_viz.png model_1.pdf
convert -limit memory 64 results/joint_qc_gmm/*/*${model_2}_bestfit_viz.png model_2.pdf
convert -limit memory 64 results/joint_qc_gmm/*/*${model_3}_bestfit_viz.png model_3.pdf
convert -limit memory 64 results/joint_qc_gmm/*/*${model_4}_bestfit_viz.png model_4.pdf
convert -limit memory 64 results/joint_qc_gmm/*/*${model_5}_bestfit_viz.png model_5.pdf
```
