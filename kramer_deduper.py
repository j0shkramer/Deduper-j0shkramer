#!/usr/bin/env python

# Documentation for a SAM file
# https://samtools.github.io/hts-specs/SAMv1.pdf

import argparse
import re


def get_args():
    parser = argparse.ArgumentParser(description="A program to take a SamTools sorted single-end SAM file and remove PCR duplicates")
    parser.add_argument("-f", "--file", help="Input SamTools-sorted single-end SAM file")
    parser.add_argument("-o", "--outfile", help="Name of the output deduplicated SAM file")
    parser.add_argument("-u", "--umi", help="File that contains all known UMIs in the SAM file, no header")
    return parser.parse_args()

args = get_args()

#########################
#   HELPER FUNCTIONS    #
#########################

def extractUMI(QNAME) -> str:
    """
    Extract the UMI from the QNAME
    """
    rname_list = QNAME.split(":")
    UMI: str = rname_list[-1]
    return UMI

# Columns: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
def extractIdentifers(samFileLine):
    """
    Extract the RNAME, bitwise flag, cigar string, chromosome, start position, and sequence from a SAM file line 
    """
    samfilelist = samFileLine.split("\t")
    qname = samfilelist[0]
    flag = int(samfilelist[1])
    chrom = int(samfilelist[2])
    startpos = int(samfilelist[3])
    cigar = samfilelist[5]
    seq_length = len(samfilelist[9])
    return qname, flag, chrom, startpos, cigar, seq_length

def generateUMISet(umifile):
    """
    Input the UMI file, and generate a set of all the UMIs for future checking"
    """
    umi_set: set = set()
    with open(umifile, "r") as umis:
        for line in umis:
            curr_umi = line.strip()
            umi_set.add(curr_umi)
    return umi_set

def getStrandedness(FLAG):
    """
    Input the bitwise flag and determine the strandedness of the read
    """
    if ((FLAG & 16)) == 16:
        return "-"
    else:
        return "+"
    
# Used to store all of the occurances in the cigar string where you should make changes to the start position
rev_add_set = set(["M", "S", "D", "N"])

def getFivePrimeStart(STARTPOS, CIGAR, STRANDEDNESS):
    """
    Determine the 5' start position of the read using its strandedness, CIGAR string and start position in the SAM file
    """
    five_prime_start_pos = STARTPOS
    # If there is no soft clipping and it is on the positive strand, the start position does not need to be movec
    if STRANDEDNESS == "+" and "S" not in CIGAR:
        return five_prime_start_pos
    # Saves the numbers from the cigar into a list
    nums_cigar = re.split(r"[A-Z]", CIGAR)
    if STRANDEDNESS == "+":
        characters_seen = 0
        for char in CIGAR:
            # The first character seen isn't S, which means there was no soft clipping on the 5' end
            if char.isalpha() and char != "S" and characters_seen == 0:
                characters_seen += 1
                return five_prime_start_pos
            # First character seen is S, meaning there was soft clipping
            elif char.isalpha() and char == "S" and characters_seen == 0:
                five_prime_start_pos -= int(nums_cigar[0])
                return five_prime_start_pos
    elif STRANDEDNESS == "-":
        characters_seen = 0
        for char in CIGAR:
            # If no characters have been seen yet and it is S, it would not change the 5' start position
            if char.isalpha() and characters_seen == 0:
                if char == "S":
                    characters_seen += 1
                else:
                    if char in rev_add_set:
                        five_prime_start_pos += int(nums_cigar[characters_seen])
                        characters_seen += 1
            elif char.isalpha():
                # Check if the character should affect the 5' start position
                if char in rev_add_set:
                    five_prime_start_pos += int(nums_cigar[characters_seen])
                    characters_seen += 1
                else:
                    characters_seen += 1
        return five_prime_start_pos

##################
#   MAIN BODY    #
##################

# Store all of the known umis
umi_set = generateUMISet(args.umi)

with open(args.file, "r") as oldSAM, open(args.outfile, "w") as newSAM: 
    # A set that will store seen reads, it will be cleared everytime we are on a new chromosome
    seen_reads = set()
    # holds the current chrom that we are on, needed for reseting the set above
    curr_chrom = ""
    for line in oldSAM:
        curr_line = line.strip()
        print(curr_line)
        if curr_line[0] == "@":
            # indicates that this is a header file
            newSAM.write(curr_line)
            newSAM.write("\n")
        else:
            qname, flag, chrom, startpos, cigar, seq_length = extractIdentifers(curr_line)
            # get the current read's unique molecular identifier
            umi = extractUMI(qname)
            # Ignore cases where UMI is not in known set
            if umi not in umi_set:
                continue
            else:
                # indicates we have advanced down to a new chromosome in the SAM file
                if chrom != curr_chrom:
                    curr_chrom = chrom
                    seen_reads.clear()
                # "+" when on forward strand, "-" when on reverse
                strandedness = getStrandedness(flag)
                # Get the current read's 5' start position
                five_prime_start_pos = getFivePrimeStart(startpos, cigar, strandedness)
                # UMI, Chrom, strand, 5' start position, sequence length
                curr_read = (umi, chrom, strandedness, five_prime_start_pos, seq_length)
                if curr_read in seen_reads:
                    continue
                else:
                    newSAM.write(curr_line)
                    newSAM.write("\n")
                    seen_reads.add(curr_read)