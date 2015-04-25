#! /usr/bin/env python
"""
Script to run blast on Genbank/EMBL files without having to first convert to fasta. 

Written by Bryan Wee.

Version 0.0.1 

"""

import sys, os
import argparse
import re
import subprocess
import filecmp

from Bio import SeqIO

actList = []

def main():    
    filesForBlasting = []
    for i in range(0,len(args.input)):

        if re.search('^\d+\.\.\d+$', args.input[i]):
            shortGenbank = cutGenbank(args.input[i-1], args.input[i])
            filesForBlasting.pop(-1) #removes original filename from the list
            filesForBlasting.append(shortGenbank) #appends the cut filename to the list
            
        elif determineFileType(args.input[i]):
            print args.input[i], "file type is:", determineFileType(args.input[i])
            filesForBlasting.append(args.input[i])

        else:
            error("Error: unable to recognize suffix of: " + args.input[i])
            
    doBlast(filesForBlasting)
    
    if args.act == True:
        loadACT(actList)

        
def doBlast(inputList):
    for i in range(0,len(inputList)):
        if determineFileType(inputList[i]) in {"genbank", "embl"}: #if file not fasta, convert to fasta
            convert2Fasta(inputList[i])
        else:
            continue
     

    for i in range(0,len(inputList) - 1): #pair up seqs for blast
        queryName = re.sub(r"\.(\w+$)", r"", inputList[i]) #strip suffixes from filename
        subjecName = re.sub(r"\.(\w+$)", r"", inputList[i+1])
        
        queryFile = ''
        subjecFile = ''
        
        actList.append(inputList[i])
                
        if determineFileType(inputList[i]) in {"genbank", "embl"}: #use fasta if converted
            queryFile = re.sub(r"\.\w+$", r".fa", inputList[i])
        else:
            queryFile = inputList[i]
            continue
        
        if determineFileType(inputList[i]) in {"genbank", "embl"}: #use fasta if converted
            subjecFile = re.sub(r"\.\w+$", r".fa", inputList[i+1])
        else:
            subjecFile = inputList[i+1]
            continue
        
        
    
        if args.blast == "blastn":
            print "Performing blastn..."
            
            #check if blast has been performed
            blast_out = queryName + ".vs." + subjecName + ".blastn.tab"
            actList.append(blast_out)
            if os.path.exists(blast_out):
                warning("Existing blast results detected, skipping...")
                pass
                            
            subprocess.Popen("blastn -query " + queryFile + ' -subject ' + subjecFile + " -outfmt 6 -out " + queryName + ".vs." + subjecName + ".blastn.tab " + args.flags, shell=True).wait()
            
        elif args.blast == "tblastx":
            print "Performing tblastx"
            
            blast_out = queryName + ".vs." + subjecName + ".tblastx.tab"
            actList.append(blast_out)
            if os.path.exists(blast_out):
                warning("Existing blast results detected, skipping...")
                pass
                            
            subprocess.Popen("tblastx -query " + queryFile + ' -subject ' + subjecFile + " -outfmt 6 -out " + queryName + ".vs." + subjecName + ".tblastx.tab", shell=True).wait()
            
    actList.append(inputList[-1]) #add last input file to ACT list
            
def loadACT(inputList):
        actCommand = " ".join(inputList)
        subprocess.Popen(["act " + actCommand + " &"], shell=True)
        #print actCommand

def convert2Fasta(file):
    convertedFilename = re.sub(r"\.\w+$", r".fa", file)
    SeqIO.convert(file, determineFileType(file), convertedFilename, 'fasta')

def determineFileType(file):
    if re.search("|".join(["gb", "gbk"]), file):
        return "genbank"
    elif re.search("|".join(["embl", "emb"]), file):
        return "embl"
    elif re.search("|".join(["fasta", "fa", "fna", "fas"]), file):
        return "fasta"
    else:
        return False
        
           
def cutGenbank(file, coordinates):
    print "Coordinates provided, cutting:", file, coordinates
    coordinatesSplit = coordinates.split('..') #splits start from stop
    startCut = int(coordinatesSplit[0])
    stopCut = int(coordinatesSplit[1])

    file2cut = open(file, "r") #assigning genbank file preceding coordinates to handle
    
    records = SeqIO.parse(file2cut, determineFileType(file))
    for rec in records:
        cutGenbank = rec[startCut:stopCut] #cuts genbank
        cutGenbank.annotations = rec.annotations
            
        cutFileName =  re.sub(r"\.(\w+$)", r".cut.\1", file) #adds "cut" before suffix
        cutFileName =  re.sub(r"cut", coordinates, cutFileName) #specifies cut in the name
        cutFile = open(cutFileName, 'w')
        SeqIO.write(cutGenbank, cutFile, determineFileType(file))
        print "Cut file written to:", cutFileName            
    
    return cutFileName
        
def error(message):
    sys.stderr.write("\n%s\n\n" % message)
    sys.exit(1)
    
def warning(message):
    sys.stderr.write("\n%s\n\n" % message)


            #function to convert genbank to fasta
#def convert2fasta(genbank, outfile, start, stop):
    #genbank = args.input[]
    #print genbank

    
#Default values    
outfile = 'blast_results.tab'
blastFlavour = 'blastn'
doACT = False
blastOptions = ''

lastarg = 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
        
    
    #takes in the input files
    parser.add_argument('input', nargs="+", action="store", help="Specify at least 2 input files")
    parser.add_argument("-o", "--output", action="store", help="Name of output file")
    parser.add_argument("-f", "--flags", action="store", help="Custom BLAST options. E.g. -f '-evalue 0.001'")
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="Verbose mode")
    parser.add_argument("-a", "--act", action="store_true", default=False, help="Run ACT")
    parser.add_argument("-b", "--blast", action="store", default="blastn", choices=("blastn", "tblastx"), help="Blast program to use")
    args = parser.parse_args()


    main() #run main script


