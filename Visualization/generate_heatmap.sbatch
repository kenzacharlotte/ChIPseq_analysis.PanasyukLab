#!/bin/bash
#SBATCH --partition=ipop-up
#SBATCH --account=metid21
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=generate_heatmap_%a.log
#SBATCH --array=0-7

module purge
module load deeptools/3.5.6

GTF="/shared/banks/genomes/mus_musculus/mm10/gtf/gencode.vM25.annotation.gtf"

# Liste des bigwig
INPUTS=( *.bw )

# Vérification
if [ ${#INPUTS[@]} -eq 0 ]; then
  echo "No bw files found"
  exit 1
fi

FILE="${INPUTS[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$FILE" .bw)

echo "Processing: $FILE"
echo "Path to the gtf : $GTF"
echo "Output name : $BASENAME"

#computeMatrix scale-regions -S "$FILE" -R "$GTF" \
#-a 1500 \
#-b 1500 \
#-o "$BASENAME.mat"

plotHeatmap -m "$BASENAME.mat" -o "$BASENAME.heatmap_zyMax.png" --yMax 15 --zMax 30

