#!/bin/bash
#SBATCH --partition=ipop-up
#SBATCH --account=metid21
#SBATCH --output=clean_sort_bam_%a.log
#SBATCH --cpus-per-task=2
#SBATCH --array=0-12

module purge
module load samtools/1.18
module load java-jdk/11.0.9.1

# List bam files
INPUTS=(*filtered.bam)

# Vérifier qu'il y a des fichiers
if [ ${#INPUTS[@]} -eq 0 ]; then
  echo "Aucun fichier *filtered.bam trouvé"
  exit 1
fi

# Récupérer les fichiers pour cette tâche
INPUT_NAME="${INPUTS[$SLURM_ARRAY_TASK_ID]}"
OUTPUT_PREFIX="${INPUT_NAME/_filtered.bam/}"

echo "Task ID: $SLURM_ARRAY_TASK_ID"

# 3 - Remove sam
#echo "3 - Remove sam & free space"
#rm ${OUTPUT_PREFIX}.sam
#echo "Done !"

# 4 - Sort bam
echo "4 - Sort bam"
samtools sort ${OUTPUT_PREFIX}_filtered.bam -o ${OUTPUT_PREFIX}_sorted.bam
echo "..."
echo "Done !"

# 5 - Mark duplicated and remove them
PICARD="picard.jar"        # Path to picard.jar file
echo "5 - Mark duplicates and remove them"
java -jar ${PICARD} MarkDuplicates \
--REMOVE_DUPLICATES true \
-I ${OUTPUT_PREFIX}_sorted.bam \
-O ${OUTPUT_PREFIX}_dup_rem.bam \
-M ${OUTPUT_PREFIX}_metrics.txt

echo "Done !"
