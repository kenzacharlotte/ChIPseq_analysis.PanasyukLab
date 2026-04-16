#!/bin/bash
#SBATCH --partition=ipop-up
#SBATCH --account=metid21
#SBATCH --output=resultat_%a.log
#SBATCH --cpus-per-task=2
#SBATCH --array=0-32

module purge
module load bowtie2/2.5.4

# Sécurise le glob (évite *_1.fq.gz si aucun fichier)
shopt -s nullglob

# Liste des fichiers R1
INPUTS=(*_1.fq.gz)

# Vérifier qu'il y a des fichiers
if [ ${#INPUTS[@]} -eq 0 ]; then
  echo "Aucun fichier *_1.fq.gz trouvé"
  exit 1
fi

# Récupérer les fichiers pour cette tâche
FILE_R1="${INPUTS[$SLURM_ARRAY_TASK_ID]}"
FILE_R2="${FILE_R1/_1.fq.gz/_2.fq.gz}"

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "R1: $FILE_R1"
echo "R2: $FILE_R2"

# Vérifier que R1 existe
if [ ! -f "$FILE_R1" ]; then
  echo "R1 introuvable"
  exit 1
fi

# Vérifier que R2 existe
if [ ! -f "$FILE_R2" ]; then
  echo "Paire manquante pour $FILE_R1"
  exit 1
fi

# Nom du sample
BASENAME=$(basename "$FILE_R1" _1.fq.gz)

# Index (IMPORTANT : préfixe, pas dossier)
INDEX="bowtie2_index/genome_index"

# Sécuriser le nombre de threads
THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "Threads: $THREADS"

# Lancement Bowtie2
bowtie2 \
  -x "$INDEX" \
  -1 "$FILE_R1" \
  -2 "$FILE_R2" \
  -p "$THREADS" \
  -S "${BASENAME}.sam"
