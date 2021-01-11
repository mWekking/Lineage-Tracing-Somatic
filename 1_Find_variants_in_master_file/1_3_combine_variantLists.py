import sys
variantDict = {}

def checkVariants(mergedVariantFile, variantDict):
    with open(mergedVariantFile) as f1:
        next(f1)
        for variant in f1:
            Chr, Loc, Ref, Sum, A, T, C, G = variant.rstrip().split("\t")

            if Chr not in variantDict:
                variantDict[Chr]    =   {Loc : {"Ref" : Ref,
                                                  "Sum" : int(Sum),
                                                  "A"   : int(A), 
                                                  "T"   : int(T), 
                                                  "C"   : int(C),
                                                  "G"   : int(G)}}


            elif Loc not in variantDict[Chr]: 
                variantDict[Chr][Loc]          = {"Ref" : Ref, 
                                                  "Sum" : int(Sum),
                                                  "A"   : int(A), 
                                                  "T"   : int(T), 
                                                  "C"   : int(C),
                                                  "G"   : int(G)}

            else:
                variantDict[Chr][Loc]["Sum"] += int(Sum)
                variantDict[Chr][Loc]["A"]   += int(A)
                variantDict[Chr][Loc]["T"]   += int(T)
                variantDict[Chr][Loc]["C"]   += int(C)
                variantDict[Chr][Loc]["G"]   += int(G)


def writeMergedFile(writingDict, outputFile):
    with open(outputFile, "w") as f1:
        f1.write("Chr\tloc\tref\tsum\tA\tT\tC\tG\n")
        for chromosome in writingDict:
            for locations in writingDict[chromosome]:
                f1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                                                    chromosome,
                                                                    locations,
                                                                    writingDict[chromosome][locations]["Ref"],
                                                                    writingDict[chromosome][locations]["Sum"],
                                                                    writingDict[chromosome][locations]["A"],
                                                                    writingDict[chromosome][locations]["T"],
                                                                    writingDict[chromosome][locations]["C"],
                                                                    writingDict[chromosome][locations]["G"]

                ))

#example: python3 1_3_combine_variantFiles.py variantsFiles.txt patient1_run123_variants.txt

#input
sampleListFile = sys.argv[1] #txt file with one relative file path per line
outputFileName = sys.argv[2]
sampleList = []
with open(sampleListFile) as f1:
    for lines in f1:
        sampleList.append(lines.rstrip())

print(sampleList)
for files in sampleList:
    checkVariants(files, variantDict)

writeMergedFile(variantDict, outputFileName)