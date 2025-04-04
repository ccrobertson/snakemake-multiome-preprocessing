This snakemake automates setup of configs and running preprocessing steps on a multiome library.

Before running update the following:
1. You may need to modify `make-config-rna.py` and `make-config-atac.py` to appropriately parse the fastq file names. The naming conventions differ over time and with different sequencing cores.
2. Run the pipeline. Add `-n` after `snakemake all` to perform a dry run.  

```
module load snakemake/7.32.4
snakemake all --jobname "{jobid}" \
    --cores 1 \
	--keep-going \
	--rerun-incomplete \
	--printshellcmds > snakemake.log 2>&1 &
```


When the nextflow pipelines have completed on all libraries, you can extract important output using the following. This will put primary results in a directory called "keep". Now you can delete the nextflow 'work' directory and free up a lot of disk space.
```
#first create the cleanup script -- creates a slurm script: work/multiome-atac/cleanup.slurm
bash scripts/cleanup-atac.sh
bash scripts/cleanup-rna.sh

#now run it (inspect it first to make sure paths are correct)
cd work/multiome-atac
sbatch cleanup.slurm > cleanup.submitted

cd work/multiome-rna
sbatch cleanup.slurm > cleanup.submitted
```

#now that you've pulled important output out of the work directory into the results directory you can delete the contents of the work dir
#to delete intermediate results files but preserve nextflow execution information (.command.log, .command.sh, etc.), you can use the following
```
cd work/multiome-atac/work
find . -not -path '*/\.*' -delete

cd work/multiome-rna/work
find . -not -path '*/.*' -delete
```


A note about version control. If you clone this repository to your workspace, the snRNAseq-nextflow and snATACseq-nextflow pipelines will not be remotely connected to the source github repositories. Below are the versions of the repositories that are currently included. For more updated versions or to keep track of changes in these repositories, replace them with clones directly from the source repositories.
* snRNAseq-nextflow (https://github.com/porchard/snRNAseq-NextFlow) - commit c4938e9f4a4db512cce707d586f53be610469b7c from Feb 10, 2025
* snATACseq-nextflow (https://github.com/porchard/snATACseq-NextFlow) - commit 18cf76233596a76c5cd16e3c14ae944a0ab0b736 from May 20, 2024
