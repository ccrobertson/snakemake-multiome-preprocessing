ROOT=${PWD}/work/multiome-rna
#ROOT is the directory from which you ran the snRNA-nextflow pipeline


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

echo -e "rsync -avL ${SOURCE}/cellbender ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/fastqc ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/interactive-barcode-rank-plots ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/multiqc ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/prune ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/qc ${DEST}" >> ${ROOT}/cleanup.slurm
echo -e "rsync -avL ${SOURCE}/starsolo ${DEST}" >> ${ROOT}/cleanup.slurm

echo -e "mv ${SOURCE} ${SOURCE}_orig" >> ${ROOT}/cleanup.slurm
echo -e "mv ${DEST} ${SOURCE}" >> ${ROOT}/cleanup.slurm
