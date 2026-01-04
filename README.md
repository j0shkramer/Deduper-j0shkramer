# Deduper

When sequences undergo PCR, there is a possibility that a template molecule is sequenced multiple times. This will artificially inflate the counts of certain reads when, in reality, their counts are not this high. If they are not removed, the information taken from the data is incorrect and can lead to improper conclusions. PCR duplicates can arise when two copies of the same original molecule get into different beads or primer lawns (depending on the sequencing technology). PCR duplicates can be identified by locating molecules that have the same UMI, 5' start position, strandness, and are located on the same chromosome.

## Running Deduper

python3 kramer_deduper.py -f (Input SAMtools-sorted SAM file) -u (File that contains all known UMIs) -o (Name of file output)
