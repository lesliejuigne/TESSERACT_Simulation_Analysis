#!/bin/bash
# Usage:
#   ./run_in_container.sh TASK_ID JOB_ID TASK_MAX JOB_NAME PARTITION JOB_TYPE MATERIAL ISOTOPE BIAS

TASK_ID=$1
JOB_ID=$2
TASK_MAX=$3
JOB_NAME=$4
PARTITION=$5
JOB_TYPE=$6
MATERIAL=$7
ISOTOPE=$8
BIAS=$9   # on | off

# --- Paths ---
BASE_DIR="/disk/data1/lze/ljuign/run3_g4sim_v11/v11-biasing"
BUILD_DIR="$BASE_DIR/TESSERACT-Octagon/build"
CONTAINER="/home/lze/ljuign/appbuilder/tes-base.sif"

# TESSERACT executable
TESSERACT_EXE="$BUILD_DIR/TesseractSim"

# Geant4 data (from conda Geant4 install)
G4DATA_HOST="/opt/conda/envs/TESinstall/share/Geant4/data"

# --- Macro + Output selection ---
if [ "$JOB_TYPE" == "internals" ]; then
    IN_DIR="$BUILD_DIR/macros/Internals/$MATERIAL"
    OUT_DIR="$BASE_DIR/internals/base"
    MACRO_FILE="$IN_DIR/gpsInMaterial_${MATERIAL}_${ISOTOPE}.mac"
    OUTPUT_FILE="$OUT_DIR/${MATERIAL}_${ISOTOPE}_${TASK_ID}_b${BIAS}.root"

elif [ "$JOB_TYPE" == "rock" ]; then
    if [ "$ISOTOPE" == "neutrons" ]; then
        IN_DIR="$BUILD_DIR/macros/LSMneutron"
        OUT_DIR="$BASE_DIR/rock/base"
        MACRO_FILE="$IN_DIR/LSM_Neutron_10M.mac"
        OUTPUT_FILE="$OUT_DIR/Rock_Neutrons_${TASK_ID}_b${BIAS}.root"
    else
        IN_DIR="$BUILD_DIR/macros/RockGammas"
        OUT_DIR="$BASE_DIR/rock/base"
        MACRO_FILE="$IN_DIR/Rock_Gamma_${ISOTOPE}_1B.mac"
        OUTPUT_FILE="$OUT_DIR/Rock_Gammas_${ISOTOPE}_${TASK_ID}_b${BIAS}.root"
    fi

elif [ "$JOB_TYPE" == "concrete" ]; then
    IN_DIR="$BUILD_DIR/macros/ConcreteGammas"
    OUT_DIR="$BASE_DIR/concrete/base"
    MACRO_FILE="$IN_DIR/Concrete_Gamma_${ISOTOPE}_1B.mac"
    OUTPUT_FILE="$OUT_DIR/Concrete_Gammas_${ISOTOPE}_${TASK_ID}_b${BIAS}.root"
fi

# --- Ensure output directory exists ---
mkdir -p "$OUT_DIR"

# --- Run inside container ---
/cvmfs/oasis.opensciencegrid.org/mis/singularity/current/bin/apptainer exec \
    --bind /disk \
    $CONTAINER bash -c "
        # --- Export required Geant4 environment variables ---
        export G4ENSDFSTATEDATA=$G4DATA_HOST/ENSDFSTATE3.0
        export G4LEDATA=$G4DATA_HOST/EMLOW8.6.1
        export G4NEUTRONHPDATA=$G4DATA_HOST/NDL4.7.1
        export G4RADIOACTIVEDATA=$G4DATA_HOST/RadioactiveDecay6.1.2
        export G4PIIDATA=$G4DATA_HOST/PII1.3
        export G4SAIDDATA=$G4DATA_HOST/SAIDDATA2.0
        export G4PARTICLEXSDATA=$G4DATA_HOST/PARTICLEXS4.1
        export G4INCLDATA=$G4DATA_HOST/INCL1.2
        export G4ABLADATA=$G4DATA_HOST/ABLA3.3
        export G4LEVELGAMMADATA=$G4DATA_HOST/PhotonEvaporation6.1
        export G4REALSURFACEDATA=$G4DATA_HOST/RealSurface2.2

        # --- Source TESSERACT environment ---
        source $BUILD_DIR/../install/bin/TesseractSimEnv.sh

        # --- Run simulation ---
        echo '=== Running task $TASK_ID on host $(hostname) ==='
        echo 'Macro:  $MACRO_FILE'
        echo 'Output: $OUTPUT_FILE'
        $TESSERACT_EXE -m $MACRO_FILE -b $BIAS -o $OUTPUT_FILE
        echo '=== Task $TASK_ID finished ==='
    "
