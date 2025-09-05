#!/bin/bash
#SBATCH --job-name=Concrete_GammasFilter            # Job name
#SBATCH --output=Concrete_GammasFilter_%A_%a.out    # Standard output log
#SBATCH --error=Concrete_GammasFilter_%A_%a.err     # Error log
#SBATCH --time=2-23:30:00
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3GB
#SBATCH --mail-user=leslie.juigne@physik.uzh.ch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-2                     # 3 tasks: 0=K40, 1=Th232, 2=U238

# Liste des isotopes
ISOTOPES=("K40" "Th232" "U238")
ISO=${ISOTOPES[$SLURM_ARRAY_TASK_ID]}

# Dossiers
INPUT_DIR="/disk/data1/lze/ljuign/G4v11_sim/Sim_Hybrid-v4.5.1/sim_data/base/rock"  # à remplacer
OUT_FOLDER="/disk/data1/lze/ljuign/G4v11_sim/Sim_Hybrid-v4.5.1/sim_data/filtered/rock"     # à remplacer
mkdir -p "$OUT_FOLDER"

# Macro ROOT
ROOT_MACRO="/disk/data1/lze/ljuign/TESSERACT_Simulation_Analysis/tesssa/cpp/VD_single_filtering.cc"

# Pattern des fichiers
BASE_PATTERN="$INPUT_DIR/Concrete_Gammas_${ISO}_*.root"

echo "----------------------------------------"
echo "Processing isotope: $ISO"
echo "Files matching pattern: $BASE_PATTERN"
echo "----------------------------------------"

for INPUT_FILE in $BASE_PATTERN; do
    [ -e "$INPUT_FILE" ] || continue
    echo "Filtering $INPUT_FILE ..."
    root -l -b -q "$ROOT_MACRO(\"$INPUT_FILE\",\"$OUT_FOLDER\")"
done

echo "Job for $ISO finished!"
