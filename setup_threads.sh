# Get the number of CPU cores
num_cores=$(nproc)

# Set the environment variable
export NUM_CORES=$num_cores
onee=1
python3 -m pip install --upgrade "jax[cuda12_local]"
export JULIA_NUM_THREADS=$((num_cores - onee)),1
