#!/usr/bin/env python
# encoding: utf-8
'''
lengthplot -- Plots the distribution of reads between 21M, 22M, 23M, and 24M, it is for mapping sRNA


@author:     Steven Hill
            
@copyright:  2013 Donald Danforth Plant Science Center. All rights reserved.
            
@license:    license

@contact:    shill@danforthcenter.org
@deffield    updated: Updated
'''

import sys
import os
from Bio import SeqIO

from optparse import OptionParser
import matplotlib

matplotlib.use('agg')
from matplotlib import pyplot
__all__ = []
__version__ = 0.2
__date__ = '2014-05-08'
__updated__ = '2014-05-08'

PROFILE = 0

def main(argv=None):
    '''Command line options.'''
    
    program_name = os.path.basename(sys.argv[0])
    program_version = "v0.2"
    program_build_date = "%s" % __updated__
 
    program_version_string = '%%Length Plot %s (%s)' % (program_version, program_build_date)
    #program_usage = '''usage: spam two eggs''' # optional - will be autogenerated by optparse
    program_longdesc = 'Length plot is a tool for plotting the length of sRNA between 21, 22, 23, and 24 basepairs' # optional - give further explanation about what the program does
    program_license = "Copyright 2014 Donald Danforth Plant Science Center                                            \
                Licensed under the Apache License 2.0\nhttp://www.apache.org/licenses/LICENSE-2.0"
 
    if argv is None:
        argv = sys.argv[1:]
        fastq, fasta, output, sam = (None,None,None,None)
        title = "Plot title"
        # setup option parser
        parser = OptionParser(version=program_version_string, epilog=program_longdesc, description=program_license)
        parser.add_option("-f", "--fastq", dest="fastq", help="input fastq file", metavar="FILE")
        parser.add_option("-F", "--fasta", dest="fasta", help="input fasta file (either this or --fastq or --sam must be specified", metavar="FILE")
        parser.add_option("-s", "--sam", dest="sam", help="input sam file (either this or --fastq or --fasta must be specified", metavar="FILE")
        parser.add_option("-o", "--output", dest="output", help="set output name (.pdf) [default: %default]", metavar="FILE")
        parser.add_option("-t", "--title", dest="title", help="set title name of the plot", metavar="STRING")

       
        # process options
        (opts, args) = parser.parse_args(argv)
        
        if opts.fastq:
            fastq = opts.fastq
        if opts.fasta:
            fasta = opts.fasta
        if opts.sam:
            sam = opts.sam
        if opts.output:
            output = opts.output
        if opts.title:
            title = opts.title
        if (fastq == None and fasta == None and sam == None) or output == None:
            indent = len(program_name) * " "
            print >> sys.stderr, program_name, "[options]", indent,  "for help use --help"
            return 2
        #lengths = count_lengths(fasta, fastq)
        # MAIN BODY #
        if sam != None:
            lengths = read_sam_lengths(sam)
        else:
            lengths = count_lengths(fasta, fastq)
        graph_length(lengths, output, title)
"""
Counts the length of each read and bins them into a dictionary

returns length dictionary     lengths = {21 : N , 22 : N, 23 : N, 24: N }

"""
def count_lengths(fasta, fastq):
    #setup
    if fasta == None and fastq == None:
        raise Exception("no input specified.")
    if fastq == None:
        fd =  open(fasta, "r")
        seq = SeqIO.parse(fd, "fasta")
    else:
        fd = open(fastq, "r")
        seq = SeqIO.parse(fd, "fastq")
        
    total = 0.0
    lengths = {21 : 0 , 22 : 0, 23 : 0, 24: 0 }
    for record in seq:
        try:
            lengths[len(record)] += 1
        except KeyError:
            pass
        finally:
            total += 1
            #potentially note
    for i in lengths:
        lengths[i] = lengths[i]/total
    fd.close()
    return lengths
"""
Reads lengths from a sam file (mapped reads only) instead of from a fastq file. returns the dictionary
of lengths which will be plotted.
"""
def read_sam_lengths(sam):
    """ sample read
    2289-2  16  scaffold01015   4825    37  9M1I10M *   0   0   TGTCAATCCTTACTATGTCT    *   XT:A:U  NM:i:2  X0:i:1  X1:i:0  XM:i:1  XO:i:1  XG:i:1  MD:Z:3A15
    """
    lengths = {21 : 0 , 22 : 0, 23 : 0, 24: 0 }
    inputfd = open(sam, "r")

    for line in inputfd:
        linearr = line.split()
        repeats = linearr[0].split("-")
        if len(repeats) < 2:
            repeats = 1
        else:
            repeats = int(repeats[1])
        #will throw exception on NaN
        try:
            curlen = int(linearr[5][:-1])
        except:
            curlen = len(linearr[7])
        if int(curlen) in lengths:
            pos = int(linearr[3])
            try:
                lengths[curlen] += 1 * repeats
            except KeyError:
                lengths[curlen] = 0
                lengths[curlen] += 1 * repeats
    print lengths
    inputfd.close()
    return lengths
"""
Expects dictionary     lengths = {21 : N, 22 : N, 23 : N, 24: N }

Builds bar plot based on the percentage of the total for each. Percent on the y-axis length on x-axis
"""
def graph_length(lengths, output, title="plot"):
    
    fig = pyplot.figure()
    plot = fig.add_subplot(1,1,1) #sets up borders
    y = [ lengths[val] for val in sorted(lengths) ]
    nbar = len(y)
    ind = range(nbar)
    labels =  ["21", "22", "23","24"]
    pyplot.ylim(0,1)
    plot.bar(ind, y, align="center")
    plot.set_ylabel('Percent of Reads')
    plot.set_title(title, fontstyle='italic')
    pyplot.xticks(ind, labels)
    pyplot.xlabel("Length (bp)")
    pyplot.savefig(output)

if __name__ == "__main__":
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'lengthplot_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    else:
        sys.exit(main())
