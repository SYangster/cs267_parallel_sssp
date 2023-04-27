#!/bin/bash

# Define the dataset you want to use
DATASET="CAL"

# Loop through the desired number of threads
for num_threads in 1 2 4 8 16 32 64; do
    echo "Running delta_stepping_parallel with $num_threads threads"
    
    # Set the environment variables
    export OMP_NUM_THREADS=$num_threads
    export OMP_PLACES=cores

    # Run the delta_stepping_parallel executable
    ./delta_stepping_parallel $DATASET
done
