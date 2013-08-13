#!/usr/bin/python
import math
import numpy
from sys import *
import matplotlib
from collections import deque
from matplotlib.lines import Line2D
test = False
server = True
def main():
    inp = None
    stepsize = 1
    titlename = None
    output = None
    windowsize = 1000
    minimum = 0
    maximum = maxint
    ymin = None
    ymax = None
    rev = None
    legend = False
    if test == True:
        inp = "plottingtest"
        output = "plottingtest.pdf"
        window = 50 # or 100ish
        display = False
    for i in range(len(argv)):
        if argv[i] in ["-i", "--input"]:
            inp = argv[i+1].split(',')
        if argv[i] in ["-f", "--title"]:
            titlename = argv[i+1]
        if argv[i] in ["-h", "--help"]:
            printHelp()
        if argv[i] in ["-s", "--step"]:
            stepsize = argv[i+1]
        if argv[i] in ["-o", "--output"]:
            output = argv[i+1]
        if argv[i] in ["-d", "--display"]:
            display = True
        if argv[i] in ["-w", "--window"]:
            windowsize = argv[i+1]
        if argv[i] in ["-m", "--window"]:
            minimum = argv[i+1]
        if argv[i] in ["-n", "--maximum"]:
            maximum = argv[i+1]
        if argv[i] in ["-y", "--ymin"]:
            (ymin) = int(argv[i+1]) 
        if argv[i] in ["-Y", "--ymax"]:
            (ymax) = int(argv[i+1]) 
        if argv[i] in ["-r", "--reverse"]:
            rev = argv[i+1].split(',')
        if argv[i] in ["-L", "--legend"]:
            legend = True
    if inp == None:
        printHelp()
        return
    if server:
        matplotlib.use('agg')
        display = False
    if windowsize != None and rev != None:
        MakeGraphWindowAvgFR(inp, rev, titlename, stepsize, windowsize, output, display, minimum, maximum, ymin, ymax, legend)
    elif windowsize != None:
        MakeGraphWindowAvg(inp, titlename, stepsize, windowsize, output, display, minimum, maximum, ymin, ymax)
    else:
        MakeGraph(inp, titlename, stepsize, output, display)

def printHelp():
    print "alignplot is a program designed to plot the output of aligndepth. That is, plot the alignment depth of different basepairs to it's location in the genome. \
alignplot supports both a single input or a list of inputs. If this is being run on a server with matplotlib, but no visual support (that is GTKAgg), or the danforth cluster, then make sure the setting at the top of the executable titled server is set to True. Using the display option is disabled when in server mode."
    print "-i    --input    <FILE,FILE,FILE,FILE>    Input files of the forward type to be graphed. These should be in CSV format and will be plotted in the positive space."
    print "-f    --title    <STRING>=titlename       The title you would like the graph to have."
    print "-h    --help                              Prints this message"
    print "-s    --step     <INT>=StepSize           The number of steps that will be taken when calculating each window. That is, a step size one will move the window one point forward."
    print "-o    --output   <INT>=Outputname         If the file is being saved as a pdf, use this flag. This will be the name of the file that is saved to disk"
    print "-d    --display                           The presence of this flag tells the program to display the plot. Display will be interactive."
    print "-w    --window   <INT>=windowsize         This is the window size to use for the sliding average. Default is 50."
    print "-m    --minimum  <INT>=minimum            Only display points from the reigon greater than or equal to minimum"
    print "-n    --maximum  <INT>=maximum            Only display points from the reigon less than or equal to maximum"    
    print "-y    --ymin     <INT>=ymin               Sets the lower bound on the y axis."
    print "-Y    --ymax     <INT>=ymax               Sets the upper bound on the y axis."
    print "-r    --reverse  <FILE,FILE,FILE,FILE>    Input files of the reverse type to be graphed. These should be in CSV format and will be plotted in the negative space."
    print "-L    --legend                            Tells the plotter to use the default 21 bp 22bp 23bp 24bp as the legend"
def MakeGraph(inp, titlename, numpoints, output, display):
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 8}

    matplotlib.rc('font', **font)
    numpoints = int(numpoints)
    alignments = {}
    with open(inp, "r") as data:
        for line in data:
            ar = line.strip().split(',')
            alignments[int(ar[0])] = ar[1]
    print "finished reading"
    if numpoints == None:
        numpoints = 10
    
    numpoints = len(alignments)
    stepsize = len(alignments)/numpoints
    print stepsize
    xpoints = []
    ypoints = []
    tot = 0
    print "sorting and building array"
    for i in sorted(alignments):
        tot += float(alignments[i])
        if int(i) % stepsize == 0:
            avg = tot/float(stepsize)
            ypoints.append(avg)
            xpoints.append(int(i))
            tot = 0
    print stepsize
    plot(xpoints, ypoints)
    if titlename == None:
        title("Align Depth/BP Pos Unmasked")
    else:
        title(titlename)
    ylabel("Alignment Depth", fontsize="8")
    xlabel("BP. Position")
    if display:
        show()
    else:
        savefig(output, dpi="300")

"""
expect forward and reverse csv?? both as list associated 1to1?
fA,fB,fC -r rA,rB,rC????
"""
def MakeGraphWindowAvgFR(inputf, inputr, titlename, stepsize, windowsize, output, display, minimum, maximum, ymin, ymax,legend):
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 8}

    if len(inputf) != len(inputr):
        raise Exception("Reverse list and Forward list are not of equal length")
    matplotlib.rc('font', **font)
    stepsize = int(stepsize)
    alignments = {}
    windowsize = int(windowsize)
    minimum = int(minimum)
    maximum = int(maximum)
    if isinstance(inputf, basestring):
        inputf = [inputf] 
    legendlist = []
    plotlist  = []
    colorlist = []
    for files in inputf:
        with open(files, "r") as data:
            for line in data:
                ar = line.strip().split(',')
                try:
                    alignments[int(ar[0])] = ar[1]
                except:
                    pass
        #setup our fields
        if stepsize > windowsize:
            raise Exception("Step size larger than window size.")
        depthdq = deque()
        locdq = deque()
        ypoints = []
        xpoints = []
        for i in range(max(alignments)):
            if i >= minimum and i <= maximum:
                if len(depthdq) < windowsize:
                    #fixes bug that was causing reigions unmapped to be grouped in the same window
                    try:
                        depthdq.appendleft(float(alignments[i]))
                        locdq.appendleft(int(i))
                    except KeyError:
                        depthdq.appendleft(0.0)
                        locdq.appendleft(int(i))
                else:
                    avg = sum(depthdq)/float(windowsize)
                    avgpos = sum(locdq)/windowsize
                    ypoints.append(avg)
                    xpoints.append(avgpos)
                    for i in range(stepsize):
                        depthdq.pop()
                        locdq.pop()
        #loading midfunction is to dodge the problem of not having the GTK backend
        from matplotlib.pyplot import plot, savefig,title, xlabel, ylabel,legend, axis
        plotlist.append(plot(xpoints, ypoints))
        if not legend:
            legendlist.append(files[:-2])
        colorlist.append(plotlist[-1][0].get_color())
    if legend:
        legendlist.append("21 bp")
        legendlist.append("22 bp")
        legendlist.append("23 bp")
        legendlist.append("24 bp")    
    #reverse it so we can use it as a queue..
    colorlist = colorlist[::-1] 
    #Now reverse it
    for files in inputr:
        with open(files, "r") as data:
            for line in data:
                ar = line.strip().split(',')
                try:
                    alignments[int(ar[0])] = ar[1]
                except:
                    pass
        #setup our fields
        if stepsize > windowsize:
            raise Exception("Step size larger than window size.")
        depthdq = deque()
        locdq = deque()
        ypoints = []
        xpoints = []
        for i in range(max(alignments)):
            if i >= minimum and i <= maximum:
                if len(depthdq) < windowsize:
                    try:
                        depthdq.appendleft(float(alignments[i]))
                    except:
                        depthdq.appendleft(0.0)
                    locdq.appendleft(int(i))
                else:
                    avg = sum(depthdq)/float(windowsize)
                    avgpos = sum(locdq)/windowsize
                    ypoints.append(avg * -1)
                    xpoints.append(avgpos)
                    for i in range(stepsize):
                        depthdq.pop()
                        locdq.pop()
        #loading midfunction is to dodge the problem of not having the GTK backend
        from matplotlib.pyplot import plot, savefig,title, xlabel, ylabel,legend, axis, axhline
        axhline(linewidth=4, color='y')

        plot(xpoints, ypoints, color=colorlist.pop())
        
    if titlename == None:
        title("Align Depth/BP Pos Unmasked Windowed Average")
    else:
        title(titlename + " Windowed Average")
    if maximum != maxint and ymin!= None and ymax != None:
        axis([minimum,maximum,ymin,ymax])
    else:
        pass
        if ymin != None and ymax != None:
            x1,x2,y1,y2 = axis()
            axis([x1,x2,ymin,ymax])
    ylabel("Alignment Depth (log10)")
    xlabel("BP. Position")
    legend(plotlist,legendlist, loc=2)
    if display:
        from matplotlib.pyplot import show
        show()
    if output != None:
        savefig(output)
        print output
    
def MakeGraphWindowAvg(input, titlename, stepsize, windowsize, output, display, minimum, maximum, ymin, ymax):
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 10}

    matplotlib.rc('font', **font)
    stepsize = int(stepsize)
    alignments = {}
    windowsize = int(windowsize)
    minimum = int(minimum)
    maximum = int(maximum)
    if isinstance(input, basestring):
        input = [input] 
    legendlist = []
    plotlist  = []
    for files in input:
        with open(files, "r") as data:
            for line in data:
                ar = line.strip().split(',')
                try:
                    alignments[int(ar[0])] = ar[1]
                except:
                    pass
        #setup our fields
        remsize = len(alignments) - int(windowsize)
        if stepsize > windowsize:
            raise Exception("Step size larger than window size.")
        depthdq = deque()
        locdq = deque()
        ypoints = []
        xpoints = []
        for i in sorted(alignments):
            if i >= minimum and i <= maximum:
                if len(depthdq) < windowsize:
                    depthdq.appendleft(float(alignments[i]))
                    locdq.appendleft(int(i))
                else:
                    avg = sum(depthdq)/float(windowsize)
                    avgpos = sum(locdq)/windowsize
                    ypoints.append(avg)
                    xpoints.append(avgpos)
                    for i in range(stepsize):
                        depthdq.pop()
                        locdq.pop()
        #loading midfunction is to dodge the problem of not having the GTK backend
        from matplotlib.pyplot import plot, savefig,title, xlabel, ylabel,legend, axis
        plotlist.append(plot(xpoints, ypoints))
        legendlist.append(files)
    if titlename == None:
        title("Align Depth/BP Pos Unmasked Windowed Average")
    else:
        title(titlename + " Windowed Average")
    if maximum != maxint:
        axis([minimum,maximum,ymin,ymax])
    else:
        pass
        x1,x2,y1,y2 = axis()
        axis([x1,x2,ymin,ymax])
    ylabel("Alignment Depth")
    xlabel("BP. Position")
    legend(plotlist,legendlist)
    if display:
        from matplotlib.pyplot import show
        show()
    if output != None:
        savefig(output)

if __name__ == "__main__":
    main()
