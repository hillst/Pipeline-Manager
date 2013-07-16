#!/usr/bin/python
import sys
import math
def main():
    output = None
    input = None
    min = 0
    max = sys.maxint
    length = 0
    lengths = None
    for i in range(len(sys.argv)): 
        if sys.argv[i] in ["-i", "--input"]:
            input = sys.argv[i+1]
        if sys.argv[i] in ["-o", "--output"]:
            output = sys.argv[i+1]
        if sys.argv[i] in ["-l", "--length"]:
            lengths = sys.argv[i+1].split(',')
        if sys.argv[i] in ["-m", "--minimum"]:
            min = int(sys.argv[i+1])
        if sys.argv[i] in ["-n", "--maximum"]:
            max = int(sys.argv[i+1])
        if sys.argv[i] in ["-h", "--help"]:
            printHelp()
        
    if output == None or input == None or lengths == None:
        printHelp()
        return    
    if input[0] == 0 or input[1] == 1:
        printHelp()
        return
    CalculateDepth(input, output, lengths, min, max)
    
"""
Re-designed to take a list of lengths and output them in the manner 21M_output 22M_output ... lengthM_output etc.
This allows the plotting step to address specifically files with some name output and look at all the files prepended
with NNM_ and plot those.

For forward and reverse building of depth files just use your own format to signal forward or reverse, and then plot them.
"""
def CalculateDepth(input, output, lengths, minimum, maximum):
    """ sample read
    2289-2  16  scaffold01015   4825    37  9M1I10M *   0   0   TGTCAATCCTTACTATGTCT    *   XT:A:U  NM:i:2  X0:i:1  X1:i:0  XM:i:1  XO:i:1  XG:i:1  MD:Z:3A15
    """
    depthDic = {}
    outputfd = {}
    inputfd = open(input, "r")
    for length in lengths:
        depthDic[int(length)] = {}
        outputfd[int(length)] = open( output +"_"+ length+"M", "w")
        
    for line in inputfd:
        linearr = line.split()
        repeats = linearr[0].split("-")
        if len(repeats) != 2:
            raise Exception("Could not split header")
        repeats = int(repeats[1]) 
        #will throw exception on NaN
        try:
            curlen = int(linearr[5][:-1])
        except:
            curlen = len(linearr[7])
        
        if not str(curlen) in lengths:
            #print "WARNING: Listed length does not equal input length, skipping..."
            pass
        else:
            pos = int(linearr[3])
            for i in range(curlen):
                try:
                    depthDic[curlen][i + pos] += 1 * repeats
                except KeyError:
                    depthDic[curlen][i + pos] = 0
                    depthDic[curlen][i + pos] += 1 * repeats
    for leng in depthDic:
        for i in sorted(depthDic[leng]):
            if i >= minimum and i <= maximum:
                outputfd[leng].write(str(i) + "," + str(math.log10(depthDic[leng][i]))+"\n")
    inputfd.close()
    for fd in outputfd:
        outputfd[fd].close()
    
def printHelp():
    print "This function is for parsing a sam file and counting the depth of each alignment. It will mark the number of reads at each \
            basepair and write to a file. This program is to be used for helping a graphing application (insuring it only needs to be done once)"
    print "-i    --input    <FILE>                       Input samfile"
    print "-o    --output   <FILE>                       Output files to be read by the graphing program. The forward set will be f_OUTPUT and the reverse set will be r_OUTPUT"
    print "-l    --length   <INT,INT,INT,INT>=lengths    Lengths of reads that will be processed. Each length with be put into it's own file with the output appended with lengthM. ex output_21M"
    print "-m    --minimum  <INT>=minimum                Beginning section of the genome to look at (default 0)"
    print "-n    --maximum  <INT>=maximum                Ending location of the genome that we will be looking at (default 9*10^19)"
    print "-h    --help                                  Prints this message"
    
if __name__ == "__main__":
    main()
