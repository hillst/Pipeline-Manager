#!/usr/bin/python
from subprocess import *
import time
import json
import threading
import Queue
from sys import argv, exit, stdout
from copy import deepcopy
from os import listdir
from os.path import isfile, join
from ctypes import ARRAY
  
def main():
    verbose = False
    strict = False
    for i in range(len(argv)):
        if argv[i] == "-h":
            helpMessage()
        if argv[i] == "-v":
            verbose = True
        if argv[i] == "-s":
            try:
                setting = argv[i+1]
            except:
                helpMessage()
                sys.exit(-1)
        if argv[i] == "-w":
            strict = True
    pipeline = Pipeline(setting, verbose, strict)
    
def helpMessage():
    print "-h                    Displays this messsage."
    print "-v                    Runs the pipeline in verbose mode."
    print "-s <settings.json>    Settings file to be input. Must be in json format."
    print "-w                    Enables the pipeline to run in strict mode. It will exit on any error."
    print ""
    print "See README for information about the settings file."
    
    
"""
Pipeline process is the general class which holds various information about each process.

procname  - Full name of the process, Name listed at the beginning of the json block.
marker    - The text that a process will output when it is finished. To add a marker to a 
            process you should wrap it in a shell script. (echo "done" at the end)
procpath  - Full path to the executable being run.
args      - Arguments to be passed to the executable
output    - The output of the process at any given time.
foreach   - The appearance of "foreach" means to run this program for each argument in
            it's value, that is, the arguments after foreach will be the last argument 
            to the program. This is specifically included for cutadapt. 
            IE: cutadapt -b ATCGG -b ATGCC inp1 
                cutadapt -b ATCGG -b ATGCC inp2
                cutadapt -b ATCGG -b ATGCC inp3
                ... etc
"""
class PipelineProcess:
    procname = ""
    marker = "done" #depcricated 
    procpath = ""
    args = []
    output = ""
    directory = ""
    dirout = ""
    priority = -1
    functionp = None
    def __init__(self, procname, procpath, args, priority, output = None, marker="done"):
        self.procname = procname
        self.procpath = procpath
        self.args = args
        self.output = output 
        self.marker = marker
        self.priority = priority
        self.process = None
    #getters and setters
    def getProcess(self): return self.process
    
    def getPriority(self): return self.priority
    
    def setPriority(self, priority): self.priority = priority
    
    def getProcname(self): return self.procname
    
    def getProcpath(self):  return self.procpath
    
    def getOutput(self): return self.output
    
    def getArgs(self): return self.args
    
    def setArgs(self, args): self.args = args
    
    def setOutput(self, output): self.output = output
        
    def getMarker(self): return self.marker
    
    def getDirectory(self): return self.directory
    
    def setDirectory(self, directory): self.directory = directory
    
    def getDirout(self): return self.dirout
    
    def setDirout(self, dirout): self.dirout = dirout
    
    def setFunctionp(self, functionp): self.functionp = functionp
    #misc functions
    """
    Wrapper function to Popen
    output expects a file object, if it exists stdout will be written to that file instead of to a pipe
    """
    def startProcessAsync(self):
        try:
            if self.output == None:
                self.process = Popen([self.getProcpath()] + self.getArgs(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
            else:
                output = open(self.output, 'w')
                self.process = Popen([self.getProcpath()] + self.getArgs(), stdin=PIPE, stdout=output, stderr=PIPE)
        except OSError:
            bcolors.printFail("Pipeline Error: the program " + self.getProcpath() + " does not exist. Exiting.")
            exit(-1)
        return self.process
    
    """ 
    getProgress calls the functionp callback. This is an instance field which contains an object to the associated function
    I'm not sure the best way to implement this. Maybe do a rough job of guessing based on inputs and things, but otherwise
    it should probably expect some kind of output from the program, and then figure out the progress that way. Again, these
    seem based on the process implementation so it will be hard to guess. 
    """
    def getProgress(self):
        self.functionp()
    
"""
Pipeline is responsible for running associated Pipeline Processes. It should execute jobs 
in order of priority. If priorities are equal, it will run them asynchroniously.
"""
class Pipeline:
    procset = {}
    jobs = Queue.Queue()
    verbose = False
    strict = False
    """
    @settings Json file containing the process and their associated fields
    """
    def __init__(self, settings, verbose = False, strict = False):
        self.settingspath = settings
        self.verbose = verbose
        self.strict = strict
        self.loadSettings()
        for prior in sorted(self.procset):
            if self.verbose == True: print "\nbeginning group of priority:", prior, "," ,len(self.procset[prior]), "processes"
            stdout.flush()
            self.startParallel(prior)
        
    """
    Begins the parallel queue for running jobs of equal priority. Spawns one process
    for each job.
    """
    def startParallel(self, prior):
        working = []
        group = self.procset[prior]
        if group[0].getDirectory() != "":
            self.loadDirectory(group[0])
        for proc in self.procset[prior]:
            proc.startProcessAsync()
            working.append(proc)
        i = 1
        for workers in working:
            workers.getProcess().wait()
            if self.verbose == True: 
                bcolors.printGood("[" + ("=" * (i * 10/len(working)) + ">" + (10 - i * 10/len(working)) * " " )+ "] " \
                                  + workers.getProcname()+": " + str(i) + "/" + str(len(working)))
                if workers.getOutput() == None: 
                    errors = workers.getProcess().stderr.readlines()
                else:
                    errors = []
                if len(errors) > 0:
                    bcolors.printFail("Errors in program: " + workers.getProcname())
                    for line in errors: 
                        print line,
                        if self.strict == True: exit(-1)
                        else: bcolors.printWarning("Use strict mode (-w) to exit on error.")
                
            if workers.getOutput() == None:    
                with open("pipeline_stats.txt", "a") as fd:
                    fd.write(workers.getProcname() + "\n")
                    fd.writelines(workers.getProcess().stdout.readlines())
                    fd.flush()
            else:
                with open("pipeline_stats.txt", "a") as fd:
                    fd.write(workers.getProcname() + "\n")
                    fd.writelines(workers.getProcess().stderr.readlines())
                    fd.flush()
            i+=1
            
    
    """
    Creates a process for each file in a directory. This function was created specifically
    for cutadapt. Copy's the parameter "process" and simply appends the filename to the end
    of the argument list. This function should be called when the directory is populated
    """
    def loadDirectory(self, process):
        path = process.getDirectory()
        #reset our process
        try:
            onlyfiles = [ f for f in listdir(path) if isfile(join(path,f)) ]
        except OSError as e:
            bcolors.printFail("Error loading directory files: " + path)
            print e
            print "Exiting"
            exit(-1)
        self.procset[process.getPriority()] = []
        try:
            dirout = process.getDirout()
        except KeyError:
            dirout = ""    
        for files in onlyfiles:
            copyp = deepcopy(process)
            args = copyp.getArgs()
            i = 0
            if copyp.getOutput() != None:
                if "$DIRECTORY" in copyp.getOutput():
                    copyp.setOutput(copyp.getOutput().replace("$DIRECTORY", dirout + files))
            for arg in args:
                if "$DIRECTORY" in str(arg):
                    args[i] = args[i].replace("$DIRECTORY", dirout + files)
                i+=1
            copyp.setArgs(args)
            try:
                self.procset[process.getPriority()].append(copyp)
                self.procset[process.getPriority()]
            except KeyError:
                self.procset[process.getPriority()] = [] 
                self.procset[process.getPriority()].append(copyp)
                
    """
    Adds a process for each item in the foreach list. It's a direct copy of the original
    process but sets the last argument of the process to be equal to one of the items in
    the foreach list.
    """
    def loadForeach(self, process, settings):
        times = settings['foreach'].split(',')
        for items in times:
            copyp = deepcopy(process)
            args = copyp.getArgs()
            args.append(items)
            copyp.setArgs(args)
            try:
                self.procset[settings['priority']].append(copyp)
            except KeyError:
                self.procset[settings['priority']] = []
                self.procset[settings['priority']].append(copyp)
                
    """
    Reads the settings json file, creates an object for each process, and then places them
    in bins according to their priority. Looks for the foreach setting first, followed
    by the directory setting, and if it finds neither it will simply add the process.
    
    foreach and directory are mutally exclusive.
    """
    def loadSettings(self):
        fd = open(self.settingspath, "r")
        sett = fd.read()
        dec = JsonSettings().decode(sett)
        for proc, settings in dec.iteritems():
            if "output" in settings:
                process = PipelineProcess(proc, settings['path'], settings['arguments'], settings['priority'], settings["output"])
            else:
                process = PipelineProcess(proc, settings['path'], settings['arguments'], settings['priority'])
            if self.verbose == True: print "adding", process.getProcname()
            if 'foreach' in settings:
                self.loadForeach(process, settings)
            elif 'directory' in settings:
                process.setDirectory(settings['directory'])
                try:
                    process.setDirout(settings['directory out'])
                except KeyError:
                    pass
                #add our dummy process
                try:
                    self.procset[settings['priority']].append(process)
                except KeyError:
                    self.procset[settings['priority']] = []
                    self.procset[settings['priority']].append(process)
            else:
                try:
                    self.procset[settings['priority']].append(process)
                except KeyError:
                    self.procset[settings['priority']] = []
                    self.procset[settings['priority']].append(process)
        fd.close()
        
        
"""Settings interface to be implemented in case a new settings file is added"""
class Settings:
    raw_settings = ""
    def __init__(self):
        pass
    """Virtual function for all child settings"""
    def decode(self, settings):
        raise NotImplementedError
    
class JsonSettings(Settings):
    def __init__(self):
        Settings.__init__(self)
    """Json Wrapper for our settings"""
    def decode(self, settings):
        try:
            return json.JSONDecoder().decode(settings)
        except ValueError as e:
            bcolors.printFail("Error parsing the settings file:")
            print e
            print "Exiting."
            exit(-1)

"""helper class (more of an enum) for coloring terminal output"""
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    @classmethod
    def printWarning(bcolors, str):
        print bcolors.WARNING + str + bcolors.ENDC
    @classmethod
    def printFail(bcolors, str):
        print bcolors.FAIL + str + bcolors.ENDC
    @classmethod
    def printGood(bcolors, str):
        print "\r" + bcolors.OKGREEN + str + bcolors.ENDC,
        stdout.flush()
if __name__ == "__main__":
    main()
