#!/bin/bash
#
#SBATCH --array=0-25
#SBATCH --time=2-23:30:00
#SBATCH --partition long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3GB
#SBATCH --mail-user=leslie.juigne@physik.uzh.ch
#SBATCH --mail-type=BEGIN,FAIL
#SBATCH --output=/disk/data1/lze/ljuign/run3_g4sim_v11/v11-biasing/slurm_jobs/temp.out

# --- Arguments ---
JOB_TYPE=$1
MATERIAL=$2
ISOTOPE=$3
BIAS=$4

# --- Determine output directory ---
if [[ "$JOB_TYPE" == "internals" ]]; then
    LOG_DIR="/disk/data1/lze/ljuign/run3_g4sim_v11/v11-biasing/internals/out"
else
    LOG_DIR="/disk/data1/lze/ljuign/run3_g4sim_v11/v11-biasing/rock/out"
fi
mkdir -p "$LOG_DIR"

# --- Run simulation and redirect stdout/stderr to the proper log ---
LOG_FILE="${LOG_DIR}/Output-${MATERIAL}-${ISOTOPE}.out"
echo "Logging to $LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "Running $JOB_TYPE $MATERIAL $ISOTOPE with bias=$BIAS"
# Example command (replace with your container call)
echo "Would run: TesseractSim -m <macro> -b $BIAS -o <output>"

# --- Record start time (only once, at first array task) ---
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo $(date +%s) > job_${SLURM_ARRAY_JOB_ID}.start
fi

# Run wrapper script (one per array task)
srun ./run_in_container.sh \
    $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_MAX \
    $SLURM_JOB_NAME $SLURM_JOB_PARTITION \
    $JOB_TYPE $MATERIAL $ISOTOPE $BIAS

# Send summary mail only at the last task
if [ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_MAX" ]; then
    start_time=$(cat job_${SLURM_ARRAY_JOB_ID}.start)
    end_time=$(date +%s)
    runtime=$((end_time - start_time))
    runtime_hms=$(printf '%02d:%02d:%02d\n' $((runtime/3600)) $((runtime%3600/60)) $((runtime%60)))

    echo -e "Simulation: $SLURM_JOB_NAME\nJob ID: $SLURM_ARRAY_JOB_ID\nPartition: $SLURM_JOB_PARTITION\nTasks: $((SLURM_ARRAY_TASK_MAX+1))\nRuntime: $runtime_hms" \
        | mail -s "SLURM job array $SLURM_JOB_NAME finished" leslie.juigne@physik.uzh.ch

    rm -f job_${SLURM_ARRAY_JOB_ID}.start
fi
