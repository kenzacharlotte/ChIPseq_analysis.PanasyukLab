#!/bin/bash
#SBATCH --partition=ipop-up
#SBATCH --account=metid21
#SBATCH --output=big_wig_%a.log
#SBATCH --cpus-per-task=2
#SBATCH --array=0-3

module purge
module load deeptools/3.5.6

CSV="/shared/projects/metid21/H3K4me3_metadata.csv"

# Vérifier que le fichier existe
if [ ! -f "$CSV" ]; then
  echo "CSV non trouvé: $CSV"
  exit 1
fi

# Ignorer le header et lire la ligne correspondant à l'array task
LINE=$(tail -n +2 "$CSV" | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")

# Extraire colonnes (séparateur = virgule)
BAM=$(echo ${LINE} | cut -d',' -f5)
SCALE=$(echo "$LINE" | cut -d',' -f6)

# Vérifier
if [ -z "$BAM" ]; then
  echo "Pas de BAM pour task $SLURM_ARRAY_TASK_ID"
  exit 1
fi

if [ -z "$SCALE" ]; then
  echo "Pas de scale factor pour $BAM"
  exit 1
fi

# Nom de sortie
OUTPUT_PREFIX=$(basename "$BAM" .bam)

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "BAM: ${BAM}"
echo "Scale factor: $SCALE"
echo "binsize = 10 for TFs"
echo "Output prefix: $OUTPUT_PREFIX"

bamCoverage \
  --bam "$BAM" \
  --outFileName "${OUTPUT_PREFIX}.bw" \
  --scaleFactor "$SCALE" \
  --normalizeUsing None \
  --numberOfProcessors 2
