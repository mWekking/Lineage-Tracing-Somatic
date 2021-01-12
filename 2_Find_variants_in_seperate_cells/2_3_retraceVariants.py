import os
import sys
    
def readMergedVariantFile(variantFileLocation):
    List_123 = []   
    sampleNames = []
    with open(variantFileLocation, "r") as f1:
        next(f1)
        for line in f1:
            chrom, loc, ref, som, A, T, C, G = line.rstrip().split("\t")
            List_123.append("{}_{}".format(chrom, loc))

    return List_123


def findVariantsInCells(mergedVariantList, individualVariantLocation, outPutFile):
    with open(outPutFile, "w") as w1:
        w1.write("sample\t" + "\t".join(mergedVariantList) + "\n")

    print("header written")
    List_123 = set(mergedVariantList)
    sampleLocation =individualVariantLocation
    allSamples = os.listdir(sampleLocation)

    doneFiles = 0
    with open(outPutFile, "a") as w1:
        for x in allSamples: 
            sampleName = "-".join(x.split("_")[1:4])
            tempList = []
            writeList = [sampleName]
            print(doneFiles)
            with open(sampleLocation + "/" + x, "r") as r1:
                next(r1)
                for line in r1:
                    chrom, loc, ref, som, A, T, C, G = line.rstrip().split("\t")
                    tempList.append("{}_{}".format(chrom, loc))
                tempList = set(tempList)

                for positions in List_123:
                    if positions in tempList:
                        writeList.append(str(1))
                    else:
                        writeList.append(str(0))

            w1.write("\t".join(writeList) + "\n")
            doneFiles += 1
#example: python3 2_3_retraceVariants.py Patient1_variants.txt variants/ Patient1_variant_table.txt


#Input
mergedRunVariants           = sys.argv[1] #all variants found in a patient (variant file generated in 1_1)
individualCellVariantFolder = sys.argv[2] #A folder with all variant files for each cell  
finalOutputName             = sys.argv[3] #output name

mergedVariantsList = readMergedVariantFile(mergedRunVariants)
findVariantsInCells(mergedVariantsList, individualCellVariantFolder, finalOutputName)