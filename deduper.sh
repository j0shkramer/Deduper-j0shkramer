#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

samfile="SAM-files/sorted_C1_SE.sam"
umis="STL96.txt"

/usr/bin/time -v ./kramer_deduper.py -f $samfile -u $umis -o "example_dedup.sam"
