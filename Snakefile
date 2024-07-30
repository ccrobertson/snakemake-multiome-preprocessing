ROOT="/scratch/scjp_root/scjp0/ccrober/preprocessing/10863-VD"

rule all:
    input: 
        "work/multiome-rna/library-config.json",
	    "work/multiome-atac/library-config.json",
	    "work/multiome-rna/nextflow.submitted",
	    "work/multiome-atac/nextflow.submitted",

rule setup_rna:
    input:
    output:
        "work/multiome-rna/library-config.json"
    params:
        root=ROOT
    shell:
	    """
	    ROOT={params.root}
	    WORK={params.root}/work/multiome-rna
	    PIPELINE={params.root}/pipelines/snRNAseq-NextFlow
	    mkdir -p ${{WORK}}
	    python scripts/make_config_rna.py ${{ROOT}}/data/fastq-rna/*.fastq.gz > ${{WORK}}/library-config.json
	    cat scripts/launch.slurm > ${{WORK}}/launch.slurm
	    echo "exec nextflow run -resume --chemistry multiome --barcode-whitelist ${{PIPELINE}}/737K-arc-v1.txt -params-file ${{WORK}}/library-config.json --results ${{WORK}}/results ${{PIPELINE}}/main.nf" >> ${{WORK}}/launch.slurm
	    """    

rule setup_atac:
    input:
    output:
        "work/multiome-atac/library-config.json"
    params:
        root=ROOT
    shell:
	    """
	    ROOT={params.root}
	    WORK={params.root}/work/multiome-atac
	    PIPELINE={params.root}/pipelines/snATACseq-NextFlow
	    mkdir -p ${{WORK}}
	    python scripts/make_config_atac.py ${{ROOT}}/data/fastq-atac/*.fastq.gz > ${{WORK}}/library-config.json
	    cat scripts/launch.slurm > ${{WORK}}/launch.slurm
	    echo "exec nextflow run -resume --chemistry multiome --barcode-whitelist ${{PIPELINE}}/737K-arc-v1.txt -params-file ${{WORK}}/library-config.json --results ${{WORK}}/results ${{PIPELINE}}/main.nf" >> ${{WORK}}/launch.slurm
	    """

rule run_rna:
	input:
	    "work/multiome-rna/library-config.json"
	output:
	    "work/multiome-rna/nextflow.submitted"
	shell:
	    """	    
	    cd work/multiome-rna
	    sbatch launch.slurm > nextflow.submitted
	    """

rule run_atac:
	input:
	    "work/multiome-atac/library-config.json"
	output:
	    "work/multiome-atac/nextflow.submitted"
	shell:
	    """	    
	    cd work/multiome-atac
	    sbatch launch.slurm > nextflow.submitted
	    """
