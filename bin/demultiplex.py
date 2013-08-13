#!/usr/bin/python
"""
author     - Steven Hill
email      - shill@danforthcenter.org
 
This is a script used for demultiplexing smallRNA self.readSection. Takes an input file with one barcode per line.
The script then parses a fasta file based on these barcodes and builds a new fasta file.

Options:
-f fasta file containing the raw data              *REQUIRED*
-b barcode file containing one barcode per line    *REQUIRED*
-k leave the barcodes on the head of the file
-t number of self.readSection to use when parsing
-l enables logging
-n keeps self.readSection with a N nucleotide label from the run data

"""
import sys
import os
import Queue
import threading
#import numpy
import time
def printHelp(errorMessage):
    print errorMessage
    print "This is a script used for demultiplexing smallRNA self.readSection. Takes an input \
     file with one barcode per line. The script then parses a fasta file based on these barcodes and builds \
     a new fasta file."
    print "-f   --fasta     <FILE>    fasta file containing the raw data              *REQUIRED*"
    print "-b   --barcodes  <FILE>    barcode file containing one barcode per line    *REQUIRED*"
    print "-B   --bcfasta   <FILE>    fasta file that has associated barcodes with the provided reads R1/R2 format"
    print "-k   --keep      leave the barcodes on the head of the file"
    print "-l   --logging   enables logging"
    print "-e   --experiment <STRING> Experiment name for the data being processesed. The program will create a new folder for the output   *REQUIRED*"
    print "-L   --label              Tells the program that the barcode file has labels associated with each barcode. They should come immediately before the barcode seperated by a space or a tab"
class DeMultiplexer:
    fasta = None
    inputsize = 0
    barcode = None
    keepBC = False
    readSection = None
    descriptors = {}
    logging = False
    nreads = False
    label = False
    barcodeList = None
    fuzzy = False
    bcfasta = False
    stats = {"badseq": 0, "corrected": 0, "indexed": 0}
    barcodeDictionary = {}
    todo = Queue.Queue()
    filehandles = []
    fastaLines = []
    numReads = None
    experiment = None
    barcodeLengths = {}
    errorList = []
    lock = threading.Lock()
    def __init__(self):
        if len(sys.argv) < 2:
            printHelp("Error: Not enough arguments.")
            sys.exit(-1)
        for i in range(len(sys.argv)):
            cur = sys.argv[i]
            if cur == "-f" or cur == "--fasta":
                self.fasta = sys.argv[i+1]
            if cur == "-b" or cur == "--barcodes":
                self.barcode = sys.argv[i+1]
            if cur == "-k" or cur == "--keep":
                self.keepBC = True
            if cur == "-t" or cur == "--threads":
                self.readSection = int(sys.argv[i+1])
                self.readSection = 1
            if cur == "-l" or cur == "--logging":
                self.logging = True
            if cur == "-B" or cur == "--bcfasta":
                self.bcfasta = True
                self.barcodefasta = str(sys.argv[i+1])
            if cur == "-e" or cur == "--experiment":
                self.experiment = str(sys.argv[i+1])
            if cur == "-L" or cur == "--label":
                self.label = True
        if self.fasta == None or self.barcode == None:
            printHelp("Fasta file or Barcode file not provided.")
            sys.exit(-1)
        if self.experiment == None:
            printHelp("Experiment name not provided. Please provide this so we can create your new directory.")
            sys.exit(-1)
        try:
            os.mkdir(self.experiment)
        except:
            pass
        statinfo = os.stat(self.fasta)
        self.inputsize = statinfo.st_size
        if self.bcfasta == True:
            self.ReadFastaWithR2()
        else:
            self.ReadFasta()

    def ReadFasta(self):
        with open(self.barcode, 'r') as fd:
            self.barcodeList = fd.readlines()
            self.labels = {}
            i = 0
            for barcode in self.barcodeList:
                if self.label:
                    spt = barcode.split()
                    if len(spt) < 2:
                        print >> sys.stderr, "--label specified, but no labels in the barcodes file, disabling"
                        self.label = False
                        break
                    self.labels[spt[1]] = spt[0].strip()
                    barcode = spt[1]
                self.barcodeList[i] = barcode.strip()
                i += 1
                self.barcodeDictionary[barcode.strip()] = []
                self.barcodeLengths[len(barcode.strip())] = True
        #open our filehandles
        self.openDescriptors() 
        #perform some number of lines at a time
        with open(self.fasta, 'r') as fd:
            i = 0
            group = []
            barcode = []
            for lines in fd:
                group.append(lines)
                if (len(group) == 4):
                    self.processReads(group, barcode)
                    group = []
                #make this optional
                comp = fd.tell()/float(self.inputsize) * 100
                if comp > i:
                    print >> sys.stderr,"\r" + ("[" + ("=" * int(comp/10)) +">"+ (" " * (10 - int(comp/10))) + "]"
                           + str(int(comp)) +  "% Complete"),
                    sys.stdout.flush()
                    i+=1
            #probably handles edge cases strangely
        #cleans up
        print >> sys.stderr, "\r" + " " * 40,
        print >> sys.stderr, "\r" +  "[==========>]100% Complete"
        self.printStats()
        self.closeDescriptors()
    """
    from start to finish handles the binning of reads based on a separate file
    with associative barcodes. To be binned, the headers must match.
    """
    def ReadFastaWithR2(self):
        with open(self.barcode, 'r') as fd:
            self.barcodeList = fd.readlines()
            i = 0
            self.labels = {}
            for barcode in self.barcodeList:
                if self.label:
                    spt = barcode.split()
                    if len(spt) < 2:
                        self.label = False
                        print >> sys.stderr, "--label specified, but no labels in the barcode file. Disabling."
                    self.labels[spt[1]] = spt[0].strip()
                    barcode = spt[1]
                self.barcodeList[i] = barcode.strip()
                i += 1
                self.barcodeDictionary[barcode.strip()] = []
        #open our filehandles
        self.openDescriptors()
        with open(self.fasta, 'r') as fd:
            with open(self.barcodefasta, 'r') as fdbc:
                count = 0
                readlines = []
                bclines = []
                mmbclen = 0
                i= 0
                linenum = 0
                for lines in fd:
                    readlines.append(lines)
                    bclines.append(fdbc.readline())
                    try:
                        if(len(readlines) == 4):
                            if bclines[0][:38] != readlines[0][:38] and bclines[0][38] \
                             == 2 and readlines[0][38] == 1:
                                #self.errorList.append(readlines)
                                self.stats["badseq"] += 1
                                barcode = None
                            else:
                                barcode = bclines[1].strip()
                                if self.keepBC == False:
                                    readlines[1] = readlines[1][len(barcode):]
                                    readlines[3] = readlines[3][len(barcode):]
                                comp = fd.tell()/float(self.inputsize) * 100
                                if comp > i:
                                    sys.stderr.write("\r" + ("[" + ("=" * int(comp/10)) +">"+ (" " * (10 - int(comp/10))) + "]"
                                        + str(int(comp)) +  "% Complete " + str(fd.tell())))
                                    sys.stderr.flush()
                                    i+=1
                                if len(barcode) != len(self.barcodeList[0]):
                                    mmbclen += 1
                                try:
                                    self.barcodeDictionary[barcode]
                                    self.stats["indexed"] += 1
                                except KeyError:
                                    barcode = None
                                    self.stats["badseq"] += 1
                                finally:
                                    self.writeBarcode(barcode, readlines)
                                    bclines = []
                                    readlines = []
                                count = 0
                            count +=1
                        linenum+=1
                    except Exception as e:
                        print >> sys.stderr, "Error found in line", linenum
                        print >> sys.stderr, e
                        readlines = []
                        bclines = []
                        count = 0
        #cleans up
        print >> sys.stderr, "\r" + " " * 40,
        print >> sys.stderr, "\r" +  "[==========>]100% Complete"
        self.printStats()
        self.closeDescriptors()     
        
    """
    Handles the processesing of each group of lines
    """
    def processReads(self,group):
        self.numReads = len(self.fastaLines)/4
        if len(group) % 4 != 0:
            print "Error, something went wrong handling the groups"
            sys.exit(-1)
        self.sortBarcodes(group)
    
    def ThreadHandler(self):
        for i in range(self.readSection):
            #calculate section...
            t = threading.Thread(target = self.FindBarcodes)
            t.daemon = True
            t.start()
        self.todo.join() #block until finished
        for barcode, lines in self.barcodeDictionary.iteritems():
            with open(self.experiment + "/" + barcode +".fq", 'wa') as fd:
                final = []
                for indicies in self.barcodeDictionary[barcode]:
                    fd.writelines(self.fastaLines[indicies:indicies+4])
                self.barcodeDictionary[barcode] = []
        with open(self.experiment + "/"+"noindex.fq", 'wa') as fd:
            for indicies in self.errorList:
                fd.writelines(self.fastaLines[indicies:indicies+4])
        #move to thread handlers, each get own file
    """
    Takes a set of 4 lines to pull the barcode from. That set is called group.
    
    If we have another file containing the associated barcodes, read it, assert the
    headers are equal, chop the barcode, then attach it to the associated list.
    """   
    def sortBarcodes(self, group):
        indexed = 0
        badseq = 0
        corrected = 0
        for length, v in self.barcodeLengths.iteritems():
            try:
                barcode = group[1][:length]
                self.barcodeDictionary[barcode] 
                group[0] = group[0].strip() + " " + "orig_bc="+barcode+" new_bc=" + barcode + " bc_diffs=0\n"
                if self.keepBC == False:
                    group[1] = group[1][length:]
                    group[3] = group[3][length:]
                indexed += 1 
            except KeyError:
                if self.nreads:
                    result = self.FixBarcodeBitwise(barcode) 
                    if result == False:
                        barcode = None
                        badseq += 1
                    else:
                        group[1] =  group[1].strip() + result
                        corrected += 1
                else:
                    #maybe append an index
                    barcode = None
                    badseq += 1
            #write the line now !!!!
            #push to writer?, barcode can also be error ;)
            self.writeBarcode(barcode, group)
            
        self.lock.acquire()
        self.stats["indexed"] += indexed
        self.stats["badseq"] += badseq
        self.stats["corrected"] += corrected
        self.lock.release()
                
    def writeBarcode(self, barcode, group):
        if barcode == None:
            fd = self.descriptors["noindex"]
        else:
            fd = self.descriptors[barcode]
        fd.writelines(group)

    def closeDescriptors(self):
        for desc, key in self.descriptors.iteritems():
            key.close()

    def openDescriptors(self):
        for barcode in self.barcodeList:
            if self.label:
                label = self.labels[barcode]
            else:
                label = barcode
            self.descriptors[barcode] = open(self.experiment + "/" + self.experiment +"_bc"+ label +".fq", 'w')
        self.descriptors["noindex"] = open(self.experiment + "/"+ self.experiment  +"noindex.fq", 'wa')
         
    """
    Checks a given broken barcode against our list and attempts
    to find the best match. It then replaces it with that match
    and sets the header accordingly
     
    Still a work in progress, not functional atm

    returns String orig_bc=n new_bc=n bc_diffs=n
    """
    def FixBarcodeBitwise(self,seq):
        distances = []
        seq_bits = self.NucToBit(seq)
        for barcode in self.barcodeList:
            barcode = barcode.strip()
            barcode_bits = self.NucToBit(barcode)
            distances.append(self.FuzzDistance(seq_bits, barcode_bits))
        if min(distances) <= 1:
            best_hit = self.barcodeList[distances.index(min(distances))]
            if seq == best_hit:
                pass
        else:
            return False
        return " orig_bc="+str(seq)+" new_bc="+str(best_hit)+" bc_diffs=" + str(min(distances)) + "\n"

    """
    We may want to use bitwise operations for a fuzzy search..
    converts given string of nucleotides to bitsi
    Finally build a numpy array (vector) and return
    """
    def NucToBit(self,orig):
        Bitwise = { "A":"11", "C":"00", "T":"10", "G":"01" }
        bitseq = ""
        for char in orig:
            try:
                bitseq += Bitwise[char]
            except:
                bitseq += Bitwise["C"]
        bitarr = numpy.array(map(int, bitseq))
        return bitarr

    """
    Finds the fuzzy distance between two given bits 
    """
    def FuzzDistance(self, a,b):
        diffs = (a != b)
        return diffs.sum()
    
    def printStats(self):
        with open("stats.txt", "w") as fd:
            total = float(sum(self.stats.values()))  
            print self.stats["badseq"]/total
            fd.write("Samples indexed: " + str(self.stats["indexed"]) + "\t"+ str(self.stats["indexed"]/total)+"\n")
            if (self.nreads):
                fd.write("Samples corrected: " + str(self.stats["corrected"]) +"\t" + str(self.stats["corrected"]/total) + "\n")
            fd.write("Bad Sequences: " + str(self.stats["badseq"]) + "\t" + str(self.stats["badseq"]/total) + "\n")

start = time.time()
demult = DeMultiplexer()
print "total time:", time.time() - start
