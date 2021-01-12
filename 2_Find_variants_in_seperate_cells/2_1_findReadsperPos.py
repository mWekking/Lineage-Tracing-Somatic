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

bamFilesList      = sys.argv[1]
refrenceFile = sys.argv[2]

with open(bamFilesList, "r") as fileList:
    for bamFile in fileList:
        bamFile = bamFile.rstrip()
        mpileupOutput  = "{}.vcf".format(bamFile.split(".")[0])

        print("Running mpileup")
        outputFinal = mpileup(bamFile, refrenceFile, mpileupOutput)

        outputFile  =  "{}_counts.txt".format(bamFile.split(".")[0])
        print("\nParsing the results")
        mpileupParser(mpileupOutput, outputFile)