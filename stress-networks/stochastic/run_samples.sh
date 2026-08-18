#!/bin/bash

for i in {1..30}; do
    echo "Running sample $i"
    /opt/homebrew/bin/python3.13 stochastic_networks.py "$i"
done


