#!/bin/bash
#SBATCH --partition=ipop-up
#SBATCH --account=metid21
#SBATCH --output=resultat_%a.log
#SBATCH --cpus-per-task=2
#SBATCH --array=0-31

module purge
module load java-jdk/11.0.9.1

# List bam files
INPUTS=(*_sorted.bam)

# Vérifier qu'il y a des fichiers
if [ ${#INPUTS[@]} -eq 0 ]; then
  echo "Aucun fichier *bam trouvé"
  exit 1
fi

# Récupérer les fichiers pour cette tâche
INPUT_NAME="${INPUTS[$SLURM_ARRAY_TASK_ID]}"
BASENAME="${INPUT_NAME/_sorted.bam/}"

echo " 6 - Mark duplicates in sample ${INPUT_NAME}"
echo "..."

java -jar picard.jar MarkDuplicates \
--REMOVE_DUPLICATES true \
-I ${BASENAME}_sorted.bam \
-O ${BASENAME}_marked_duplicates.bam \
-M ${BASENAME}_metrics.txt 

echo "Done !"
