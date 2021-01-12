import subprocess
import sys
import os

def mpileup(bamFile, refrenceFile, outputFile):
    subprocess.call(["samtools", 
                    "mpileup", 
                    "-A",
                    "-f",  refrenceFile,
                    bamFile,
                    "-o", outputFile,
                    "--ff", "DUP"])
    
    return outputFile

def mpileupParser(outputFile, countFile):
    subprocess.call(["perl", 
                    "pileup2base.pl", 
                    outputFile, 
                    "0", 
                    countFile])

#example usage: python3 2_1_findReadsperPos.py bamList.txt hg19_withchr.fa

#input
bamFilesList      = sys.argv[1] #a list of all bam files, one file per line
refrenceFile      = sys.argv[2] ##hg19.fa chromosomes names should be chr1,chr2 etc.. not only a digit

with open(bamFilesList, "r") as fileList:
    for bamFile in fileList:
        bamFile = bamFile.rstrip()
        mpileupOutput  = "{}.vcf".format(bamFile.split(".")[0])

        print("Running mpileup")
        outputFinal = mpileup(bamFile, refrenceFile, mpileupOutput)

        outputFile  =  "{}_counts.txt".format(bamFile.split(".")[0])
        print("\nParsing the results")
        mpileupParser(mpileupOutput, outputFile)