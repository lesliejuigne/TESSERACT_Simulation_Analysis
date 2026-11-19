# STABILITY CHECK SCRIPT FOR TESSERACT SIMULATIONS
# =================================================
# This script runs TesseractSim simulations using the base TESSERACT geometry. 
# It performs stability checks by running simulations with identical and different seeds,
# comparing the outputs to ensure the version installed is stable.

# To run this script: 
# 1. Ensure you have TesseractSim built and the path set correctly.
# 2. Source the G4 environment.
# 3. Execute the script using the command : python3 SimTester.py

# =================================================

# Modules
import subprocess
import os
import uproot
import numpy as np

# Define paths to the executable and macros
path = "path/to/the/build/folder"
executable = os.path.join(path, "TesseractSim")

# =================================================
class PyRunner:
    def __init__(self, file_path):
        self.file_path = file_path  

    def get_root_tree(self):
        with uproot.open(self.file_path) as root_file:
            tree = root_file["events"]
            Params = list(tree.keys())
            rx   = tree[Params[7]].array(library="np")
            ry   = tree[Params[8]].array(library="np")
            rz   = tree[Params[9]].array(library="np")
            t    = tree[Params[13]].array(library="np")
            edep = tree[Params[16]].array(library="np")
        return rx, ry, rz, t, edep

def run_macro(executable, macro_file, output_file):
    cmd = [executable, "-m", macro_file, "-o", output_file]
    print(f"\nRunning: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"Finished: {output_file}")

def summarize_data(file_path):
    rx, ry, rz, t, edep = PyRunner(file_path).get_root_tree()
    return {
        "n_events": len(edep),
        "edep_mean": np.mean(edep),
        "edep_std": np.std(edep),
        "x_mean": np.mean(rx),
        "y_mean": np.mean(ry),
        "z_mean": np.mean(rz),
    }

def compare_runs(f1, f2, expect_identical=True):
    s1 = summarize_data(f1)
    s2 = summarize_data(f2)

    print(f"\nComparing {f1} vs {f2}:")
    all_ok = True

    for key in s1.keys():
        same = np.isclose(s1[key], s2[key], rtol=1e-9, atol=1e-9)
        if expect_identical and not same:
            print(f"  {key}: DIFFERENT (expected same) — {s1[key]:.4e} vs {s2[key]:.4e}")
            all_ok = False
        elif not expect_identical and same:
            print(f"  {key}: SAME (expected different) — {s1[key]:.4e}")
            all_ok = False
        else:
            print(f"  {key}: OK ({s1[key]:.4e} vs {s2[key]:.4e})")

    return all_ok


def delete_files(files):
    for f in files:
        if os.path.exists(f):
            os.remove(f)
            print(f"Deleted {f}")
        else:
            print(f"{f} not found, skipping")

# =================================================

if __name__ == "__main__":
    runs = [
        {"macro": os.path.join(path, "macros/StabilityCheck/setup_test_1.mac"), "output": "test_1.0.root"},
        {"macro": os.path.join(path, "macros/StabilityCheck/setup_test_1.mac"), "output": "test_1.1.root"},
        {"macro": os.path.join(path, "macros/StabilityCheck/setup_test_2.mac"), "output": "test_2.0.root"},
    ]

    for run in runs:
        if not os.path.exists(run["output"]):
            run_macro(executable, run["macro"], run["output"])
        else:
            print(f"Skipping {run['output']} (already exists)")

    # ============== Checks ==============
    print("\n=== Test 1: identical seeds (1.0 vs 1.1) ===")
    ok1 = compare_runs("test_1.0.root", "test_1.1.root", expect_identical=True)
    print("Result:", "PASS" if ok1 else "FAIL")

    print("\n=== Test 2: different seeds (1.0 vs 2.0) ===")
    ok2 = compare_runs("test_1.0.root", "test_2.0.root", expect_identical=False)
    print("Result:", "PASS" if ok2 else "FAIL")

    # ============== Cleanup ==============
    delete_after = input("\nDelete generated ROOT files? (y/n): ").strip().lower()
    if delete_after == "y":
        delete_files([r["output"] for r in runs])
        print("All test files deleted — simulations verified successfully. Have fun with the simulations!")
    else:
        print("Files kept, something might be wrong.")
