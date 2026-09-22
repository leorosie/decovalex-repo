#!/bin/bash

CASE_DIR="cases_yaml"

for sample in $(seq 2 30); do
    for yaml_file in "$CASE_DIR"/*.yaml; do
        echo "Running sample $sample with $yaml_file"
        python3.13 stochastic_networks_with_stress.py "$sample" "$yaml_file"
    done
done