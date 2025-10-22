#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log


samfile="/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam"
umis="/projects/bgmp/joshkram/bioinfo/Bi624/Deduper-j0shkramer"

/usr/bin/time -v python3 kramer_deduper.py -f $samfile -u $umis -o "example_dedup.sam"
