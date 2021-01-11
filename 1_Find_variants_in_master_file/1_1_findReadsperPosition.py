import subprocess
import sys
import os

def mpileup(bamFile, refrenceFile, outputFile):
    subprocess.call(["samtools", 
                    "mpileup",
                    "-A",  
                    "-f",  refrenceFile,
                    bamFile,
                    "-o", outputFile]) #, "--ff", "DUP"
    
    return outputFile

def mpileupParser(outputFile, countFile):
    subprocess.call(["perl", 
                    "pileup2base.pl", 
                    outputFile, 
                    "0", 
                    countFile])

#example: python3 1_1_findReadsperPosition.py sorted_run_1_cell_1.bam hg19_withchr.fa run1_cell1_counts.txt

#Input files 
bamFile         = sys.argv[1]
refrenceFile    = sys.argv[2] #hg19.fa chromosomes names should be chr1,chr2 etc.. not only a digit
parsedCountFile = sys.argv[3] #countFileName

#mpileup output name
mpileupOutput  =  "{}.vcf".format(bamFile.split(".")[0])

#running mpileup and parsing the data
print("Running mpileup")
outputFinal = mpileup(bamFile, refrenceFile, mpileupOutput)
print("\nParsing the results")
mpileupParser(mpileupOutput, parsedCountFile)
