#!/bin/bash
#SBATCH --job-name=SSiFilter
#SBATCH --output=SSiFilter_%A_%a.out
#SBATCH --error=SSiFilter_%A_%a.err
#SBATCH --time=2-23:30:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=long
#SBATCH --mem-per-cpu=3GB
#SBATCH --mail-user leslie.juigne@physik.uzh.ch
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --array=0-11   # 12 tasks: 0=Ac228, 1=Bi214, 2=Co60, 3=Cs137, 4=K40, 5=Pb210, 6=Pb212, 7=Pb214, 8=Th232, 9=Tl208, 10=U235, 11=U238


ISOTOPES=("Ac228" "Bi214" "Co60" "Cs137" "K40" "Pb210" "Pb212" "Pb214" "Th232" "Tl208" "U235" "U238")
ISO=${ISOTOPES[$SLURM_ARRAY_TASK_ID]}

INPUT_DIR="/disk/data1/lze/ljuign/G4v11_sim/Sim_Hybrid-v4.5.1/sim_data/base/internals"
OUT_FOLDER="/disk/data1/lze/ljuign/G4v11_sim/Sim_Hybrid-v4.5.1/sim_data/filtered/internals"

# Macro ROOT
ROOT_MACRO="/disk/data1/lze/ljuign/TESSERACT_Simulation_Analysis/tesssa/cpp/VD_single_filtering.cc"

# Pattern pour tous les fichiers de l'isotope
BASE_PATTERN="$INPUT_DIR/SSi_${ISO}_*.root"

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
