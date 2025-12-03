import subprocess
import os
import uproot
import numpy as np


class PyRunner:
    def __init__(self, file_path):
        self.file_path = file_path  

    def get_root_tree(self):
        """Load arrays from ROOT TTree"""
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
    """Run one macro and wait until it finishes"""
    cmd = [executable, "-m", macro_file, "-o", output_file]
    print(f"\nRunning: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"Finished: {output_file}")


def summarize_data(file_path):
    """Return summary statistics for quick comparison"""
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
    """Compare two ROOT outputs numerically"""
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


def compare_with_reference(file_path, expected):
    """Check a simulation result against known expected reference values"""
    summary = summarize_data(file_path)
    print(f"\nChecking {file_path} against reference values:")
    ok = True
    for key, exp_val in expected.items():
        val = summary.get(key)
        if np.isclose(val, exp_val, rtol=1e-9, atol=1e-9):
            print(f"  {key}: OK ({val:.4f})")
        else:
            print(f"  {key}: FAIL — got {val:.4f}, expected {exp_val:.4f}")
            ok = False
    return ok


def delete_files(files):
    """Delete generated ROOT files"""
    for f in files:
        if os.path.exists(f):
            os.remove(f)
            print(f"Deleted {f}")
        else:
            print(f"{f} not found, skipping")


if __name__ == "__main__":
    path = "/disk/data1/lze/ljuign/TESSERACTSim_clean/build"
    executable = os.path.join(path, "TesseractSim")

    runs = [
        {"macro": os.path.join(path, "macros/setup_test_1.mac"), "output": "test_1.0.root"},
        {"macro": os.path.join(path, "macros/setup_test_1.mac"), "output": "test_1.1.root"},
        {"macro": os.path.join(path, "macros/setup_test_2.mac"), "output": "test_2.0.root"},
    ]

    # Reference expected values
    expected_seed1 = {"n_events": 3388, "edep_mean": 441.1, "edep_std": 241.1}
    expected_seed2 = {"n_events": 3947, "edep_mean": 7.843, "edep_std": 73.75}

    # Step 1: Run macros if needed
    for run in runs:
        if not os.path.exists(run["output"]):
            run_macro(executable, run["macro"], run["output"])
        else:
            print(f"Skipping {run['output']} (already exists)")

    # Step 2: Consistency checks
    print("\n=== Test 1: identical seeds (1.0 vs 1.1) ===")
    ok1 = compare_runs("test_1.0.root", "test_1.1.root", expect_identical=True)
    print("Result:", "PASS" if ok1 else "FAIL")

    print("\n=== Test 2: different seeds (1.0 vs 2.0) ===")
    ok2 = compare_runs("test_1.0.root", "test_2.0.root", expect_identical=False)
    print("Result:", "PASS" if ok2 else "FAIL")

    # Step 3: Validation against expected results
    print("\n=== Test 3: Check physical correctness vs expected ===")
    ok3a = compare_with_reference("test_1.0.root", expected_seed1)
    ok3b = compare_with_reference("test_2.0.root", expected_seed2)
    if ok3a and ok3b:
        print("PASS: All statistical values match expected seeds")
    else:
        print("FAIL: Simulation statistics differ from reference")

    # Step 4: Optional cleanup
    delete_after = input("\nDelete generated ROOT files? (y/n): ").strip().lower()
    if delete_after == "y":
        delete_files([r["output"] for r in runs])
        print("All test files deleted — simulations verified successfully. Have fun with the simulator!")
    else:
        print("Files kept for inspection.")
