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
    
    #outputFile = "H005cl07ChIP1_refrence.vcf"
    #outputFile  = bamFile[:-4] + "_refrence.vcf"
    return outputFile

def mpileupParser(outputFile, countFile, inDelFile):
    subprocess.call(["perl", 
                    "pileup2base.pl", 
                    outputFile, 
                    "0", 
                    countFile])




bamFile      = sys.argv[1]
refrenceFile = sys.argv[2]
outputFile  =  "{}_refrence.vcf".format(bamFile.split("_")[0])
countFile  = sys.argv[3]

print("Running mpileup")
outputFinal = mpileup(bamFile, refrenceFile, outputFile)

inDelFile  = "inDels"
print("\nParsing the results")
mpileupParser(outputFile, countFile, inDelFile)
