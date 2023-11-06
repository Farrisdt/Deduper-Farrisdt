#!/usr/bin/env python
import re
import argparse

#argparse
def get_args():
    parser = argparse.ArgumentParser(description="Takes in a sorted Sam file and returns a Sam file with all repeat reads removed. The first encounter of each read is sent to the deduped output file. Subsequent reads are sent to the duplicates file. Unrecognized UMIs are sent to an error file. Also reports number of unique reads per chromosome in a report file.")
    parser.add_argument("-f", "--file", help="Path to sorted Sam file", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="Name of output file", type=str, required=True)
    parser.add_argument("-u", "--umi", help="Path to file holding list of UMIs", type=str, required=True)
    return parser.parse_args()
args = get_args()

#Variables
infile = args.file
outfile = args.outfile
dupfile = f"{outfile[:-4]}_duplicates.txt" #Name for duplicates file
errorfile = f"{outfile[:-4]}_error.txt" #Name for file holding rejected reads
reportfile = f"{outfile[:-4]}_report.txt" #Name for file holding rejected reads
umifile = args.umi
dict = {} #used to hold reads we have seen, per chromosome
umidict = {} #holds uniqe umis
currchrom = "0" #current chromosome, char due to X and Y
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
x=1 #used to track number of reads for report.
with open(infile, "r") as file, open(outfile, "w") as out, open(dupfile, "w") as dup, open(errorfile, "w") as error, open(reportfile, "w") as report:
    for lines in file:
        if lines[0]=="@":
            #write out headers
            out.write(lines)
        else:
            line = lines.split("\t")
            qname = line[0].split(":")
            if not currchrom == line[2]:
                #clears dict for each new chromosome to save space
                dict.clear()
                if not currchrom=="0":
                    r = currchrom+":"+str(x)+"\n"
                    report.write(r)
                currchrom = line[2]
                if first: #if the first line, the currchrom is set to default.
                    first=False #will break for loop in the future, avoid dropping chrom 1.
                else:
                    x=0
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
            key = umi + str(truepos) + strand
            if key in dict.keys():
                dup.write(lines)
            else:
                out.write(lines)
                dict[key] = 0
                x=x+1
    r = currchrom+":"+str(x)+"\n"
    report.write(r)
