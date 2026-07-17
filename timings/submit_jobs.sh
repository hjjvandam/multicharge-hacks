#!/bin/bash
for run in {1..10}; do
  sbatch dftd4-2000.sh
done
for run in {1..10}; do
  sbatch dftd4-4000.sh
done
for run in {1..10}; do
  sbatch dftd4-8000.sh
done
