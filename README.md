This snakemake automates setup of configs and running preprocessing steps on a multiome library.

Before running update the following:
1. Set the ROOT path at top of the Snakefile.
2. You may need to modify `make-config-rna.py` and `make-config-atac.py` to appropriately parse the fastq file names. The naming conventions differ over time and with different sequencing cores.
3. Run the pipeline. Add `-n` after `snakemake all` to perform a dry run.  

```
module load snakemake/7.32.4
snakemake all --jobname "{jobid}" \
    --cores 1 \
	--keep-going \
	--rerun-incomplete \
	--printshellcmds > snakemake.log 2>&1 &
```


When the snATAC-nextflow pipeline has completed on all libraries, you can extract important output using the following. This will put primary results in a directory called "keep". Now you can delete the nextflow 'work' directory and free up a lot of disk space.
```
#first create the cleanup script -- creates a slurm script: work/multiome-atac/cleanup.slurm
bash scripts/cleanup-atac.sh

#now run it (inspect it first to make sure paths are correct)
cd work/multiome-atac
sbatch cleanup.slurm
```
