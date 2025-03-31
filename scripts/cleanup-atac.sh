ROOT=${PWD}/work/multiome-atac
#ROOT is the directory from which you ran the snATAC-nextflow pipeline
#Example: ROOT=/scratch/scjp_root/scjp0/ccrober/preprocessing/7799-VD/work/multiome-atac

SOURCE=${ROOT}/results
DEST=${ROOT}/keep

echo -e "#!/bin/bash" > ${ROOT}/cleanup.slurm
echo -e "#SBATCH --time=16:00:00" >> ${ROOT}/cleanup.slurm
echo -e "#SBATCH --mem=1000M" >> ${ROOT}/cleanup.slurm
echo -e "#SBATCH --account=scjp1" >> ${ROOT}/cleanup.slurm
echo -e "#SBATCH --output=%u-%x-%j.log" >> ${ROOT}/cleanup.slurm
echo -e "#SBATCH --error=%u-%x-%j.err" >> ${ROOT}/cleanup.slurm
echo -e "#SBATCH --mail-user=ccrober@umich.edu" >> ${ROOT}/cleanup.slurm
echo -e "#SBATCH --mail-type=END,FAIL" >> ${ROOT}/cleanup.slurm
echo -e "#SBATCH --signal=B:TERM@60" >> ${ROOT}/cleanup.slurm
echo -e "mkdir -p ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/ataqv ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/bigwig ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/counts ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/macs2 ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/mark_duplicates ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/multiqc ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/prune ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/plot-barcodes-matching-whitelist ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/fragment-file ${DEST}" >> ${ROOT}/cleanup.slurm

echo -e "mv ${SOURCE} ${SOURCE}_orig" >> ${ROOT}/cleanup.slurm
echo -e "mv ${DEST} ${SOURCE}" >> ${ROOT}/cleanup.slurm
