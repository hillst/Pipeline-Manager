#!/usr/bin/python
import sys
from sys import argv
from os import listdir
from os.path import isfile, join, exists
"""
Helper script that expects all the files to plot to be in their own directory. Takes each CSV in this directory
and draws them on a single plot.
"""
def printHelp():
    print "Helper script that expects all files to plot to exist in their own directory. Takes each CSV in this directory and draws them on a single plot"
    print ""
    print "-n   --maximum   <INT>   Maximum position on the genome to plot"
    print "-m   --minimum   <INT>   Minimum position on the genome to plot"

samples = "barcodes"
srcfiles = "plots"
base = "./alignplot -s 20 -w 100 "
bclist = []
toexec = []
minimum = None
maximum = None
for i in range(len(argv)):
    if argv[i] == "-n" or argv[i] == "--maximum":
        maximum = argv[i+1]
    if argv[i] == "-m" or argv[i] == "--minimum":
        minimum = argv[i+1]
    if argv[i] == "-h" or argv[i] == "--help":
        printHelp()    
if minimum != None:
    base += "-m " + minimum + " "
if maximum != None:
    base += "-n " + maximum + " "
base += "-i "

with open(samples,'r') as fd:
    for bc in fd: bclist.append(bc.strip())
for bc in bclist:
    path=srcfiles+"/"+bc
    final = ""
    onlyfiles = [ f for f in listdir(path) if isfile(join(path,f)) ]
    forward = []
    reverse = []
    for file in onlyfiles:
        if file[-2:] == "_f":
            forward.append(file)
        if file[-2:] == "_r":
            reverse.append(file)
    #removing last comma each time
    for file in forward: final += path+"/"+file +","
    final = final[:-1]
    final += " -r "
    for file in reverse: final += path+"/"+file +","
    final = final[:-1]
    toexec.append(base + final + " -f " + bc + " -o" + " " + bc + "_plots_21m24m_test.pdf")
for todo in  toexec: print todo
