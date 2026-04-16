#!/bin/bash
#SBATCH --partition=ipop-up
#SBATCH --account=metid21
#SBATCH --output=resultat_%a.log
#SBATCH --cpus-per-task=2
#SBATCH --array=0-12

module purge
module load bedtools/2.31.1

BLACKLIST="mm10.blacklist_ENCODE2016.bed"

# List bam files
INPUTS=(*_dup_rem.bam)

# Vérifier qu'il y a des fichiers
if [ ${#INPUTS[@]} -eq 0 ]; then
  echo "Aucun fichier *bam trouvé"
  exit 1
fi

# Récupérer les fichiers pour cette tâche
INPUT_NAME="${INPUTS[$SLURM_ARRAY_TASK_ID]}"
BASENAME="${INPUT_NAME/_dup_rem.bam/}"

echo " 8 - Exclude blacklist regions from ${BASENAME}_dup_rem.bam"
echo "..."
bedtools intersect -v -abam ${BASENAME}_dup_rem.bam \
-b ${BLACKLIST} > ${BASENAME}_blacklistremoved.bam
echo "Done !"
