ROOT=$1
#ROOT is the directory from which you ran the snATAC-nextflow pipeline
#Example: ROOT=/scratch/scjp_root/scjp0/ccrober/preprocessing/7799-VD/work/multiome-atac

SOURCE=${ROOT}/results
DEST=${ROOT}/keep


echo -e "mkdir -p ${DEST}" > ${ROOT}/launch-cleanup.sh
echo -e "rsync -avL ${SOURCE}/ataqv ${DEST}" >> ${ROOT}/launch-cleanup.sh
echo -e "rsync -avL ${SOURCE}/bigwig ${DEST}" >> ${ROOT}/launch-cleanup.sh
echo -e "rsync -avL ${SOURCE}/macs2 ${DEST}" >> ${ROOT}/launch-cleanup.sh
echo -e "rsync -avL ${SOURCE}/multiqc ${DEST}" >> ${ROOT}/launch-cleanup.sh
echo -e "rsync -avL ${SOURCE}/prune ${DEST}" >> ${ROOT}/launch-cleanup.sh
echo -e "rsync -avL ${SOURCE}/plot-barcodes-matching-whitelist ${DEST}" >> ${ROOT}/launch-cleanup.sh
echo -e "rsync -avL ${SOURCE}/fragment-file ${DEST}" >> ${ROOT}/launch-cleanup.sh
