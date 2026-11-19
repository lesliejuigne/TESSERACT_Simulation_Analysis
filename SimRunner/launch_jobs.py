import os
import subprocess
import job_config
# Path to the SLURM submit script
submit_script = "./submit_job.sh"


# Loop over internals
try:
    from job_config import internals, bias_options
    for material, isotopes in internals.items():
        for iso in isotopes:
            for bias in bias_options:
                cmd = f"sbatch --job-name='{material}_{iso}_{bias}' {submit_script} internals {material} {iso} {bias}"
                print(f"Submitting: {cmd}")
                subprocess.run(cmd, shell=True)
except Exception as e:
    print(f"Error processing internals: {e}")
    
# Loop over rock gammas
try:
    from job_config import rock_gammas,bias_options
    for iso in rock_gammas:
        for bias in bias_options:
            cmd = f"sbatch --job-name='rock_{iso}_{bias}' {submit_script} rock placeholder {iso} {bias}"
            print(f"Submitting: {cmd}")
            subprocess.run(cmd, shell=True)
except Exception as e:
    print(f"Error processing rock gammas: {e}")

# Loop over concrete gammas
try: 
    from job_config import concrete_gammas,bias_options
    for iso in concrete_gammas:
        for bias in bias_options:
            cmd = f"sbatch --job-name='concrete_{iso}_{bias}' {submit_script} concrete placeholder {iso} {bias}"
            print(f"Submitting: {cmd}")
            subprocess.run(cmd, shell=True)
except Exception as e:
    print(f"Error processing concrete gammas: {e}")
    
# Rock neutrons
try: 
    from job_config import rock_neutrons,bias_options
    if rock_neutrons:
        for bias in bias_options:
            cmd = f"sbatch --job-name='rock_neutrons_{bias}' {submit_script} rock placeholder neutrons {bias}"
            print(f"Submitting: {cmd}")
            subprocess.run(cmd, shell=True)
except Exception as e:
    print(f"Error processing rock neutrons: {e}")