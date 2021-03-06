#!/usr/bin/env python
# encoding: utf-8
'''
pcorrelation -- Plots the regression coefficient lane comparing the depth of alignment at each base pair from one sample to another.

For each bp, compare the alignment depth from one sample to another ( two samples aligned to the same reference ).

@author:     Steven Hill
            
@copyright:  2013 Donald Danforth Plant Science Center. All rights reserved.
            
@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import sys
import os
import math

from optparse import OptionParser
import matplotlib

matplotlib.use('agg')
from matplotlib import pyplot
__all__ = []
__version__ = 0.1
__date__ = '2013-08-13'
__updated__ = '2013-08-13'

PROFILE = 0

def main(argv=None):
    '''Command line options.'''
    
    program_name = os.path.basename(sys.argv[0])
    program_version = "v0.1"
    program_build_date = "%s" % __updated__
 
    program_version_string = '%%Pearson Correlation Plot %s (%s)' % (program_version, program_build_date)
    #program_usage = '''usage: spam two eggs''' # optional - will be autogenerated by optparse
    program_longdesc = '''''' # optional - give further explanation about what the program does
    program_license = "Copyright 2013 Donald Danforth Plant Science Center                                            \
                Licensed under the Apache License 2.0\nhttp://www.apache.org/licenses/LICENSE-2.0"
 
    if argv is None:
        argv = sys.argv[1:]
        inputx, inputy, output = (None,None,None)
        title = "Plot title"
        # setup option parser
        parser = OptionParser(version=program_version_string, epilog=program_longdesc, description=program_license)
        parser.add_option("-x", "--inputx", dest="inputx", help="Input file for the first sample, expects a csv file where column 0 is the position, column 1 is the depth.", metavar="FILE")
        parser.add_option("-y", "--inputy", dest="inputy", help="Input file for the second sample, expects a csv file where column 0 is the position, column 1 is the depth.", metavar="FILE")
        parser.add_option("-u", "--unifiedinput", dest="uinput", help="Unified input expects one input file that is either comma, tab, or space delemited. Expects 2 or 3 columns. If there are 2 columns, it will use column 0 as sample 1 and column 1 as sample 2. If there are 3 columns, it will use column 1 and column 2 instead.")
        parser.add_option("-o", "--output", dest="output", help="set output name (.pdf) [default: %default]", metavar="FILE")
        parser.add_option("-t", "--title", dest="title", help="set title name of the plot", metavar="STRING")

        
        # process options
        (opts, args) = parser.parse_args(argv)
        
        if opts.inputx:
            inputx = opts.inputx
        if opts.inputy:
            inputy = opts.inputy
        if opts.output:
            output = opts.output
        if opts.title:
            title = opts.title
        if opts.uinput:
            uinput = opts.uinput
        if ((inputy == None or inputx == None) and uinput == None)or output == None:
            indent = len(program_name) * " "
            print >> sys.stderr, program_name, "[options]", indent,  "for help use --help"
            return 2
        #lengths = count_lengths(fasta, fastq)
        # MAIN BODY #
        if uinput == None:
            newinputs = read_inputs(inputx, inputy)
        else:
            newinputs = read_inputs(-1, -1, uinput=uinput)
        rvalue = pearson_def(newinputs[0], newinputs[1])
        if uinput == None: 
            graph_correlation("Sample1", "Sample2" ,newinputs[0], newinputs[1], rvalue, output, title)
        else:
            graph_correlation(inputx, inputy, newinputs[0], newinputs[1], rvalue, output, title)

"""
Accepts 2 filenames as input, expects them to be csv coordinate files

builds 2 lists where the index is the genomic position and the value is the depth at that position
returns a list for each as a tuple (inputx, inputy)

"""
def read_inputs(inputx, inputy, uinput = None):
    if uinput is None:
        xpos, ypos = {}, {}
        with open(inputx, 'r') as fd:
            for line in fd:
                linespl = line.strip().split(',')
                xpos[int(linespl[0])] = int(linespl[1])
        with open(inputy, 'r') as fd:
            for line in fd:
                linespl = line.strip().split(',')
                ypos[int(linespl[0])] = int(linespl[1])
        maximum = max(max(ypos.keys()), max(xpos.keys())) #boom tuple
        fxlist = range(maximum)
        fylist = range(maximum)
        for i in xrange(maximum):
            try:
                fxlist[i] = xpos[i]
            except KeyError:
                fxlist[i] = 0.0
            try:
                fylist[i] = ypos[i]
            except KeyError:
                fylist[i] = 0.0
    else:
        i = 0
        fxlist, fylist = [],[]
        with open(uinput, 'r') as fd:
            for line in fd:
                linesp = line.split(",")
                if len(linesp) == 2 : # tab/space
                    raise Exception("invalid file format, should be csv")
                if len(linesp) == 3: #bogus first
                    try:
                        fxlist.append(float(linesp[1]))
                        fylist.append(float(linesp[2]))
                    except:
                        print >> sys.stderr, line
                elif(len(linesp) == 2):
                    try:
                        fxlist.append(float(linesp[0]))
                        fylist.append(float(linesp[1]))
                    except:
                        print >> sys.stderr, line
    return fxlist, fylist

def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

def pearson_def(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in xrange(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)
"""
expects two lists, returns slope and intercept (slope, intercept)
  Step 3: Find ΣX, ΣY, ΣXY, ΣX2.
            ΣX = 311 
            ΣY = 18.6 
            ΣXY = 1159.7 
            ΣX2 = 19359 

  Step 4: Substitute in the above slope formula given.
            Slope(b) = (NΣXY - (ΣX)(ΣY)) / (NΣX2 - (ΣX)2)
            = ((5)*(1159.7)-(311)*(18.6))/((5)*(19359)-(311)2)
            = (5798.5 - 5784.6)/(96795 - 96721)
            = 13.9/74
            = 0.19 

  Step 5: Now, again substitute in the above intercept formula given.
            Intercept(a) = (ΣY - b(ΣX)) / N 
            = (18.6 - 0.19(311))/5
            = (18.6 - 59.09)/5
            = -40.49/5
            = -8.098
"""
def regression_def(x, y):
    if len(x) != len(y):
        raise Exception("x and y not of the same length")
    sumx = sum(x)
    sumy = sum(y)
    sumxy = sum([ x[i] * y[i] for i in range(len(x)) ])
    sumxsq = sum( [ val * val for val in x ] )
    #cast to float to avoid possible integer division
    slope = (len(x) * sumxy - sumx*sumy)/ float(len(x) * sumxsq - sumx * 2)
    intercept = (sumy - slope * sumx) / float(len(x))
    return slope,intercept


"""
Expects dictionary     lengths = {21 : N, 22 : N, 23 : N, 24: N }

Builds bar plot based on the percentage of the total for each. Percent on the y-axis length on x-axis
"""
def graph_correlation(xname, yname, xlist, ylist, rvalue, output, title="plot"):
    
    fig = pyplot.figure()
    plot = fig.add_subplot(1,1,1) #sets up borders
    
    #labels =  ["21M", "22M", "23M","24M"]
    plot.scatter(xlist, ylist)
    plot.set_ylabel(yname)
    plot.set_title(title + " r=" + str(rvalue), fontstyle='italic')
    pyplot.xlabel(xname)
    yregline = [ rvalue * val for val in xrange(len(xlist)) ]
    xregline = xrange(len(xlist))
    #plot.plot(xregline, yregline)
    maximum = max(int(max(xlist)), int(max(ylist)))
    regress = regression_def(xlist, ylist)
    yreg = [ i * regress[0] + regress[1] for i in range(int(max(xlist))) ]
    xreg = [i for i in range(int(max(xlist)))]
    x1,x2,y1,y2 = plot.axis()
    plot.plot( xreg, yreg, color="r" )
    plot.axis((0,x2,0,y2))
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
