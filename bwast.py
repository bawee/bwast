#! /usr/bin/env python
"""
Script to run blast on Genbank/EMBL files without having to first convert to fasta. 

Written by Bryan Wee.

Version 0.0.1 - 20150425

"""

import sys, os
import argparse
import re
import subprocess
from argparse import RawTextHelpFormatter

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
            if args.verbose: print args.input[i], "file type is:", determineFileType(args.input[i])
            filesForBlasting.append(args.input[i])
            


        else:
            error("Error: unable to recognize suffix of: " + args.input[i])
            
    doBlast(filesForBlasting)
    
    if args.act == True:
        loadACT(actList)

        
def doBlast(inputList):
    for file in inputList:
        if determineFileType(file) in {"genbank", "embl"}: #if file is genbank, convert to fasta
            convert2Fasta(file)
        else:
            pass
        
    
    for i in range(0,len(inputList) - 1): #pair up seqs for blast
        queryName = re.sub(r"\.(\w+$)", r"", inputList[i]) #strip suffixes from filename
        subjecName = re.sub(r"\.(\w+$)", r"", inputList[i+1])
        
        queryName = os.path.basename(queryName)
        subjecName = os.path.basename(subjecName)
        
        
        queryFile = ''
        subjecFile = ''
        
        actList.append(inputList[i])
                
        if determineFileType(inputList[i]) in {"genbank", "embl"}: #Query sequence for BLAST. if genbank, use fasta file generated earlier
            queryFile = re.sub(r"\.\w+$", r".fa", inputList[i])
        else:
            #queryFile = inputList[i]
            
            #check if query is is a multi fasta
            filetype = determineFileType(inputList[i])
            records = list(SeqIO.parse(inputList[i], filetype))
            if len(records) > 1:
                if args.verbose: print "Multifasta detected, concatenating prior to BLASTing"
                mergedFile = mergeRecords(inputList[i])
                queryFile = mergedFile
                    
            else:
                queryFile = inputList[i]
                pass
        
        if determineFileType(inputList[i+1]) in {"genbank", "embl"}: #Subject sequence for BLAST, if genbank, use fasta file generated earlier
            subjecFile = re.sub(r"\.\w+$", r".fa", inputList[i+1])
        else:
            
            #check if query is is a multi fasta 
            filetype = determineFileType(inputList[i+1])
            records = list(SeqIO.parse(inputList[i+1], filetype))
            if len(records) > 1:
                if args.verbose: print "Multifasta detected, concatenating prior to BLASTing"
                mergedFile = mergeRecords(inputList[i+1])
                subjecFile = mergedFile
                    
            else:
                subjecFile = inputList[i+1]
                pass
        
        blastType = args.blast
        
        if args.verbose: print "Performing blast"
        blastOptionsPre = (args.flags if args.flags else "")
        blastOptions = re.sub(r"-(\w+)\s", r"\1_", blastOptionsPre)
        blastOptions = re.sub(r"\s+", r".", blastOptions)
        if args.verbose: print "with options: %s" % (blastOptionsPre)
        
        blast_out = "%s.vs.%s.%s.%s.tab" % (queryName, subjecName, blastOptions, blastType)
        
        actList.append(blast_out) #append blast file to ACT input list
        
        if os.path.exists(blast_out): #check if blast output exists
            if args.verbose: warning("Existing blast results detected, skipping...")
            pass
        
        #run BLAST
        subprocess.Popen("%s -query %s -subject %s -outfmt 6 -out %s %s" % (blastType, queryFile, subjecFile, blast_out, blastOptionsPre), shell=True).wait()
        #print "%s -query %s -subject %s -outfmt 6 -out %s %s" % (blastType, queryFile, subjecFile, blast_out, blastOptionsPre)


    actList.append(inputList[-1]) #add last input file to ACT list
            
def loadACT(inputList):
        actCommand = " ".join(inputList)
        subprocess.Popen(["act " + actCommand + " &"], shell=True)
        
def mergeRecords(file): #adapted from SeqHandler by NF Alikhan (github.com/happykhan/seqhandler)

#SeqHandler is a script for merging, converting and splitting sequence files (Genbank, EMBL, fasta and others). Please use it to merge multi-Genbank files before running bwast.py

    filetype = determineFileType(file) #determine file type
    readInMultifasta = open(file, "r")
    records = list(SeqIO.parse(readInMultifasta, filetype))

    mergingFile = records[0]
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    contigs = SeqFeature(FeatureLocation(0, len(mergingFile) ), type="fasta_record",\
                strand=1)
    contigs.qualifiers["note"] = records[0].name #pull out contig number of first contig
    mergingFile.features.append(contigs) #append first contig to mergingFile
    for nextRecord in records[1:]:
        contigs = SeqFeature(FeatureLocation(len(mergingFile), len(mergingFile) + len(nextRecord)), type="fasta_record",\
                strand=1)
        contigs.qualifiers["note"] = nextRecord.name 
        mergingFile.features.append(contigs) #append subsequent contigs to mergingFile
        mergingFile += nextRecord
    mergingFile.name = records[0].name
    mergingFile.description = records[0].description
    mergingFile.annotations = records[0].annotations

    for feature in mergingFile.features:
        if feature.type == 'source':
            mergingFile.features.remove(feature)
    contigs = SeqFeature(FeatureLocation(0, len(mergingFile)), type="source", strand=1)
    mergingFile.features.insert(0,contigs)
    merged_file = re.sub(r"\.\w+$", r".merged.fa", file)
    out_handle = open(merged_file, "w")
    SeqIO.write(mergingFile, out_handle, filetype)
    return merged_file

#END of code adapted from SeqHandler.

def convert2Fasta(file):
    convertedFilename = re.sub(r"\.\w+$", r".fa", file)
    SeqIO.convert(file, determineFileType(file), convertedFilename, 'fasta')

def determineFileType(file):
    if re.search("|".join(["gb$", "gbk$"]), file):
        return "genbank"
    elif re.search("|".join(["embl$", "emb$"]), file):
        return "embl"
    elif re.search("|".join(["fasta$", "fa$", "fna$", "fas$"]), file):
        return "fasta"
    else:
        return False
        
           
def cutGenbank(file, coordinates):
    if args.verbose: print "Coordinates provided, cutting:", file, coordinates
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
        if args.verbose: print "Cut file written to:", cutFileName            
    
    return cutFileName
        
def error(message):
    sys.stderr.write("\n%s\n\n" % message)
    sys.exit(1)
    
def warning(message):
    sys.stderr.write("\n%s\n\n" % message)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
Wrapper script to run blast on Genbank/EMBL files without having to first convert to fasta.
    
Input: Genbank, EMBL or Fasta
    
Requires: BLAST+, ACT and BioPython on your PATH
    ''', formatter_class=RawTextHelpFormatter)
        
    
    #takes in the input files
    parser.add_argument('input', nargs="+", action="store", help="Specify at least 2 input files")
    parser.add_argument("-f", "--flags", action="store", help="Custom BLAST options, enclosed in quotes. E.g. -f '-task blastn -evalue 0.001'")
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="Verbose mode")
    parser.add_argument("-a", "--act", action="store_true", default=False, help="Run ACT.")
    parser.add_argument("-b", "--blast", action="store", default="blastn", choices=("blastn", "tblastx"), help="Blast program to use. Default [blastn]")
    args = parser.parse_args()


    main() #run main script


