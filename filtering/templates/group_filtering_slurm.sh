#!/bin/bash
#SBATCH --job-name=CuFilter # Job name
#SBATCH --output=CuFilter_%A_%a.out # Standard output and error log
#SBATCH --error=CuFilter_%A_%a.err # Error log
#SBATCH --time=2-23:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3GB
#SBATCH --partition=long
#SBATCH --mail-user leslie.juigne@physik.uzh.ch # Replace with your email
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --array=0-2   # 3 tasks: 0=K40, 1=Th232, 2=U238 (Adjust based on the number of isotopes)


ISOTOPES=("K40" "Th232" "U238") # List of isotopes to process
ISO=${ISOTOPES[$SLURM_ARRAY_TASK_ID]}

INPUT_DIR="path/to/your/simulation/data"  # Replace with your actual input directory
OUT_FOLDER="path/to/your/output/data"  # Replace with your desired output directory

MAX_FILES=80 # Maximum number of files to process per job
ROOT_MACRO="TESSERACT_Simulation_Analysis/tesssa/cpp/VD_group_filtering.cc" # Path to your ROOT macro

BASE="Cu_${ISO}" # Base name for the dataset

echo "----------------------------------------"
echo "Processing dataset: $BASE (0–100)"
echo "----------------------------------------"

root -l -b -q "$ROOT_MACRO(\"$INPUT_DIR/$BASE\", $MAX_FILES, \"$OUT_FOLDER\")"

echo "Job for $BASE finished!"
