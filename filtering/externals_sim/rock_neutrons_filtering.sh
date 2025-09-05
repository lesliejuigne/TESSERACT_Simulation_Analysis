#!/bin/bash
#SBATCH --job-name=Rock_Neutrons_TotFilter            # Job name
#SBATCH --output=Rock_Neutrons_TotFilter.out    # Standard output log
#SBATCH --error=Rock_Neutrons_TotFilter.err     # Error log
#SBATCH --time=2-23:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3GB
#SBATCH --partition=long
#SBATCH --mail-user=leslie.juigne@physik.uzh.ch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0


# Dossiers
INPUT_DIR="/disk/data1/lze/ljuign/G4v11_sim/Sim_Hybrid-v4.5.1/sim_data/base/rock"  # à remplacer
OUT_FOLDER="/disk/data1/lze/ljuign/G4v11_sim/Sim_Hybrid-v4.5.1/sim_data/filtered/rock"     # à remplacer
mkdir -p "$OUT_FOLDER"

# Macro ROOT
ROOT_MACRO="/disk/data1/lze/ljuign/TESSERACT_Simulation_Analysis/tesssa/cpp/VD_single_filtering.cc"

# Pattern des fichiers
BASE_PATTERN="$INPUT_DIR/Rock_Neutrons_Tot_*.root"

echo "----------------------------------------"
echo "Files matching pattern: $BASE_PATTERN"
echo "----------------------------------------"

for INPUT_FILE in $BASE_PATTERN; do
    [ -e "$INPUT_FILE" ] || continue
    echo "Filtering $INPUT_FILE ..."
    root -l -b -q "$ROOT_MACRO(\"$INPUT_FILE\",\"$OUT_FOLDER\")"
done


echo "Job for Neutrons finished!"
