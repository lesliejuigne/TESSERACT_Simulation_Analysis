#!/bin/bash
#SBATCH --job-name=PCuFilter
#SBATCH --output=PCuFilter_%A_%a.out
#SBATCH --error=PCuFilter_%A_%a.err
#SBATCH --time=2-23:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3GB
#SBATCH --partition=long
#SBATCH --mail-user leslie.juigne@physik.uzh.ch
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --array=0-2   # 3 tasks: 0=K40, 1=Th232, 2=U238


ISOTOPES=("K40" "Th232" "U238")
ISO=${ISOTOPES[$SLURM_ARRAY_TASK_ID]}

INPUT_DIR="/disk/data1/lze/ljuign/G4v11_sim/Sim_Hybrid-v4.5.1/sim_data/base/internals"
OUT_FOLDER="/disk/data1/lze/ljuign/G4v11_sim/Sim_Hybrid-v4.5.1/sim_data/filtered/internals"

# Macro ROOT
ROOT_MACRO="/disk/data1/lze/ljuign/TESSERACT_Simulation_Analysis/tesssa/cpp/VD_single_filtering.cc"

# Pattern pour tous les fichiers de l'isotope
BASE_PATTERN="$INPUT_DIR/PCu_${ISO}_*.root"

echo "----------------------------------------"
echo "Processing isotope: $ISO"
echo "Files matching pattern: $BASE_PATTERN"
echo "----------------------------------------"

# Boucle sur tous les fichiers correspondant
for INPUT_FILE in $BASE_PATTERN; do
    [ -e "$INPUT_FILE" ] || continue
    echo "Filtering $INPUT_FILE ..."
    root -l -b -q "$ROOT_MACRO(\"$INPUT_FILE\",\"$OUT_FOLDER\")"
done

echo "Job for $ISO finished!"