This snakemake automates setup of configs and running preprocessing steps on a multiome library.

Before running update the following:
1. Set the ROOT path at top of Snakefile.
2. You may need to modify make-config-rna.py and make-config-atac.py to appropriately parse the fastq file names. The naming conventions differ over time and with different sequencing cores.


```
module load snakemake/7.32.4
snakemake all --jobname "{jobid}" \
    --cores 1 \
	--keep-going \
	--rerun-incomplete \
	--printshellcmds > snakemake.log 2>&1 &
```


Dry-run
```
snakemake all -n --jobname "{jobid}" \
    --cores 1 \
	--keep-going \
	--rerun-incomplete \
	--printshellcmds > snakemake.log 2>&1 &
```