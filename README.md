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
