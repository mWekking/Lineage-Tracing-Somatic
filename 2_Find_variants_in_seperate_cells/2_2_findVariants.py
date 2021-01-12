import sys
from itertools import takewhile, repeat
import subprocess
import os

def rawincount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

def readCountFile(mergedTXT, writeFile, totalLines):
    foundPositionsCount = 0
    with open(mergedTXT, "r") as r1, open(writeFile, "w") as w1:
        readHeader  = r1.readline().rstrip().split("\t")
        w1.write("chr\tloc\tref\tsum\tA\tT\tC\tG\n")
        
        for nOL, line in enumerate(r1):
            #chrom, loc, ref, A, T, C, G, a, t, c, g, insertion, deletion = line.rstrip().split("\t")
            chrom, loc, ref, A, T, C, G, a, t, c, g = line.rstrip().split("\t")

            if nOL % 1000 == 0:
                #print("working on chromosome {}".format(chrom), end = "\r")
                #print("\n")
                print("Positions Processed {} of {}".format(nOL, totalLines), end = "\r")
            Acounts = int(A) + int(a)
            Tcounts = int(T) + int(t)
            Ccounts = int(C) + int(c)
            Gcounts = int(G) + int(g)

            totalCounts = [Acounts, Tcounts, Ccounts, Gcounts]

            if totalCounts.count(0) >= 3:
                pass
            else:
                foundPositionsCount += 1
                w1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom  ,
                                                             loc    ,
                                                             ref,
                                                             sum(totalCounts),
                                                             Acounts,
                                                             Tcounts,
                                                             Ccounts,
                                                             Gcounts))

    return foundPositionsCount     

#example: python3 2_2_findVariants_sort.py 2_2_FileList.txt 

#input
inputFileList  = sys.argv[1] #a file with one count file per line

with open(inputFileList, "r") as r1:
    for countFile in r1:
        countFile  = countFile.rstrip()
        outputFile = "{}_variants.txt".format(countFile[0:-11])

        totalNumberOfLines = rawincount(countFile)
        print(totalNumberOfLines)
        numberOfVariants = readCountFile(countFile, outputFile, totalNumberOfLines)
        print("number of variants found: {}".format(numberOfVariants))
