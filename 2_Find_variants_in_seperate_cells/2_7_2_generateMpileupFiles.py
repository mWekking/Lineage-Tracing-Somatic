import subprocess
import sys
import os

def mpileup(bamFile, refrenceFile, coordinates, outputFile):
    subprocess.call(["samtools", 
                    "mpileup", 
                    "-A",
                    "-l", coordinates,
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

    os.remove(outputFile)



bamFile      = sys.argv[1]
refrenceFile = sys.argv[2]
coordinatesFile = sys.argv[3]
outputFile  = "{}_refrence.vcf".format(bamFile.split(".")[0])

print("Running mpileup")
outputFinal = mpileup(bamFile, refrenceFile, coordinatesFile, outputFile)

countFile  = bamFile.split(".")[0] + "_counts.txt" #the name of the variant output file
print("\nParsing the results")
mpileupParser(outputFile, countFile)

