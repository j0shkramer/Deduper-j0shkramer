**Define the Problem**

When sequences undergo PCR, there is a possibiliy that a template molecule is sequenced multiple times. This will artifically inflate the counts of certain reads when, in reality, their counts are not this high. If they are not removed, the information taken from the data is incorrect, and can lead to improper conclusions. PCR duplicates can arise when two copies of the same orginal molecule get into different beads or primer lawns (depends on the sequencing technology). PCR duplicates can be identifed by locating molecules that have the same UMI, 5' start position, strandness and are located on the same chromosome.

**Psuedocode**

- Go through the UMI input file and put all the UMIs into a dictonary 
- Open the orginal SAM file and create a new SAM file to copy the non-PCR duplicates to
- Go line-by-line through the orginal SAM file
- If the line begins with an "@" copy it to the new SAM file and add a new line character
    - This means it is from the header section and should be kept
- Else, it is an alignment section, and should be seperated into a list object
    - Extract the RNAME, bitwise flag, cigar string, chromosome, start position, and sequence
    - Extract the UMI from the RNAME
    - Check if the UMI is apart of the known UMIs
        - If it is not, continue to the next sequence
    - Check the bitwise flag to see if the the sequence is reverse complemented or not 
        - If it is reversed complemented we will have to change the start position to get the 5' start position
    - Check the cigar string to see if there was soft clipping
        - This is dependent on whether the sequence is reserved complemented or not
            - If it reserve complemented, we have to see if the "S" is at the start of the cigar string versus at the start otherwise
        - Store the number in front of the "S"
        - if S in cigar string?
    - Determine the 5' start position if there was soft clipping or if it was reverse complemented
        - If it was on the forward strand and there was reverse complementing, subtract the start position by soft clipping
        - If it is on the reverse strand, add the length of the sequence and the soft clipping to the start position
    - Create a list with the format [UMI, 5' start position, chromosome, Reverse or forward strand]
    - Check if this list is already in a dictonary 
        - If it is not, add copy the line to the new SAM file, add a new line character, and create a new key in the dictonary
        - If it is already in the dictonary, continue 

**Functions**

The purpose of most of these functions is to make the main function more readible 

```
findSoftClipping(str cigarString, bool isReverseComplement) -> int:
    """
    If you have determined that a sequence has soft clipping, this function will extract the number of bases that it has been soft-clipped by
    """
    Input: 10S50M5S, True
    Output: 5
```

```
extractUMI(str RNAME) -> str:
    """
    Extract the UMI from the RNAME
    """
    Input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:AACGCCAT
    Output: AACGCCAT
```

```
extractIdentifers(str samFileLine) -> str, str, str, str, int, str
    """
    Extract the RNAME, bitwise flag, cigar string, chromosome, start position, and sequence from a SAM file line 
    """
    Input: NS:AAC	0	1	15	36	4M	*	0	0	TCCA 6AEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
    Output: NS:AAC, 0, 4M, 1, 15, TCCA
```