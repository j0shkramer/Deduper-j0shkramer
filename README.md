# Deduper

PCR amplification can cause a single template molecule to be sequenced multiple times. This leads to **PCR duplicates**, which artificially inflate read counts and can bias downstream analyses.

PCR duplicates typically occur when multiple copies of the same original molecule enter different beads or primer lawns during sequencing. These duplicates can be identified because they share all of the following:

- **Unique Molecular Identifier (UMI)**
- **5′ start position**
- **Strand**
- **Chromosome**

Removing duplicates ensures that read counts reflect true biological abundance rather than PCR artifacts.

---

## Requirements

- `python ≥ 3.10` environment

**SAM File**

- Should be previously sorted using Samtools

**File Containing UMIs**

- Should be a text file, with each known UMI listed on its own line

---

## Running Deduper

```bash
python3 kramer_deduper.py \
    -f <input SAMtools-sorted SAM file> \
    -u <file containing all known UMIs> \
    -o <output file name>
