#!/bin/bash
#SBATCH --partition=ipop-up
#SBATCH --account=metid21
#SBATCH --output=peak_calling.macs2_%a.log
#SBATCH --cpus-per-task=2
#SBATCH --array=0-4

module purge
module load macs2/2.2.7.1

# Input commun
INPUT=final_bam/n2_I3_final_sorted.bam

# Liste des IP
IP=($(find final_bam -name "*n2_*[PV][1-4]*sorted.bam" ! -name "*n2_P4_*"))

# Vérification
if [ ${#IP[@]} -eq 0 ]; then
  echo "Aucun fichier trouvé"
  exit 1
fi

# Sélection par job
IP_FILE=${IP[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "$IP_FILE" .bam)

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "IP: $IP_FILE"
echo "Input: $INPUT"

# === BROAD ===
echo "Calling broad peaks..."

macs2 callpeak \
  -t "$IP_FILE" \
  -c "$INPUT" \
  -f BAMPE \
  -g mm \
  --pvalue 0.03 \
  -n "$BASENAME" \
  --broad \
  --nomodel \
  --extsize 200 \
  --outdir "macs2_n2_12h_${BASENAME}_broad"

echo "Done broad"

# === NARROW ===
echo "Calling narrow peaks..."

macs2 callpeak \
  -t "$IP_FILE" \
  -c "$INPUT" \
  -f BAMPE \
  -g mm \
  --pvalue 0.03 \
  -n "$BASENAME" \
  --nomodel \
  --extsize 200 \
  --outdir "macs2_n2_12h_${BASENAME}_narrow"

echo "Done narrow"
