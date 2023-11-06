#!/usr/bin/env python
import re
import argparse

#argparse
def get_args():
    parser = argparse.ArgumentParser(description="Takes in a sorted Sam file and returns a Sam file with all repeat reads removed. The first encounter of each read is sent to the deduped output file. Subsequent reads are sent to the duplicates file. Unrecognized UMIs are sent to an error file. Also reports number of unique reads per chromosome in a report file.")
    parser.add_argument("-f", "--file", help="Path to sorted Sam file", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="Name of output file", type=str, required=True) #not needed here
    parser.add_argument("-u", "--umi", help="Path to file holding list of UMIs", type=str, required=True)
    return parser.parse_args()
args = get_args()

#Variables
infile = args.file
outfile = args.outfile
errorfile = f"{outfile[:-4]}_error.txt" #Name for file holding rejected reads
tempfile = "./dedupertempfile.sam"
umifile = args.umi
umidict = {} #holds uniqe umis
first = True #used to avoid error on reading in first chromosome

#functions
def isReverse(flag) -> bool:
    '''Checks to see if string is forward or reverse oriented. Returns true if reversed.'''
    if ((int(flag) & 16) == 16):
        return True #read is reverse
    else:
        return False

def findLength(cigar) -> int:
    '''Returns length of read including deletions and skipped regions. Does not include insertions or clipping.'''
    length = 0

    #parse for desired sections
    match = re.findall(r'(\d+)M', cigar)
    delete = re.findall(r'(\d+)D', cigar)
    skip = re.findall(r'(\d+)N', cigar)
    lengthList = match + delete + skip

    #add sections together, cast as int
    for item in lengthList:
        length = length + int(item)
    return length

def truePosition(cigar, pos, flag) -> int:
    '''Returns the corrected position for a read, taking into account soft clipping and reverse orientation.'''
    newpos = 0
    fwrd = 0 #lefthand soft clipping
    rvse = 0 #righthand soft clipping
    #parse cigar for clipping
    snip = re.findall(r'(\d+)S', cigar)
    #assign snip numbers based on which side they are on. Accounts for one sided snipping.
    if len(snip)>0:
        if len(snip) == 2:
            fwrd = snip[0]
            rvse = snip[1]
        elif cigar[len(cigar)-1]=="S":
            rvse = snip[0]
        else:
            fwrd = snip[0]
    # + or - strand check
    if not isReverse(flag):
        # if + strand subtract first clip from pos
        newpos = int(pos) - int(fwrd)
    else:
        #if - strand add length of section and snipping
        length = findLength(cigar)
        newpos = int(pos) + int(length) + int(rvse)
        newpos-=1 #remove 1 due to 0 index error
    return newpos

#code
with open(umifile, "r") as file:
    for umi in file:
        umi = umi.strip()
        umidict[umi] = 0

with open(infile, "r") as file, open(tempfile, "w") as out, open(errorfile, "w") as error:
    for lines in file:
        if lines[0]=="@":
            #write out headers
            out.write(lines)
        else:
            lines=lines.strip()
            line = lines.split("\t")
            qname = line[0].split(":")
            flag = line[1]
            pos = line[3]
            cigar = line[5]
            umi = qname[len(qname)-1]
            strand = "x"
            if not umi in umidict.keys():
                error.write(lines)
                continue
            if isReverse(flag):
                strand = "-"
            else:
                strand = "+"
            truepos = truePosition(cigar, pos, flag)
            qname[3] = truepos
            newqname = ":".join([str(z) for z in qname])
            line[0] = newqname
            oldpos = f"{pos}\n"
            line.append(oldpos)
            newline = "\t".join([str(z) for z in line])
            out.write(newline)
