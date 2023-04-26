#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET=$1

./delta_stepping_parallel $DATASET
./dijkstra $DATASET
./delta_stepping_serial $DATASET
