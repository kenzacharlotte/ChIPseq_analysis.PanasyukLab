#!/bin/bash
#SBATCH --partition=ipop-up
#SBATCH --account=metid21
#SBATCH --output=resultat_%a.log
#SBATCH --cpus-per-task=2
#SBATCH --array=0-32

module purge
module load samtools/1.18

# List bam files
INPUTS=(*.sam)

# Vérifier qu'il y a des fichiers
if [ ${#INPUTS[@]} -eq 0 ]; then
  echo "Aucun fichier *sam trouvé"
  exit 1
fi

# Récupérer les fichiers pour cette tâche
INPUT_NAME="${INPUTS[$SLURM_ARRAY_TASK_ID]}"
OUTPUT_PREFIX="${INPUT_NAME/.sam/}"

echo "Task ID: $SLURM_ARRAY_TASK_ID"

# 1 - Convert SAM to BAM
echo "1 - Start converting : $INPUT_NAME"
echo "Output prefix will be : $OUTPUT_PREFIX"
echo "..."

echo "Convert SAM to BAM" 
samtools view -h -b ${OUTPUT_PREFIX}.sam > ${OUTPUT_PREFIX}.bam
echo "Done !"

# 2 - Filter reads with mapping quality higher than 40
echo "2 - Filter reads with mapping quality higher than 40"
echo "Output prefix will be : $OUTPUT_PREFIX_filtered"
samtools view -b ${OUTPUT_PREFIX}.bam -q 40 > ${OUTPUT_PREFIX}_filtered.bam
echo "Done !"

