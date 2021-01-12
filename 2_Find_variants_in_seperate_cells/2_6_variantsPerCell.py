import os
import sys
def getSNPs(SNPFile):
    variationCount = 0
    SNPcount       = 0
    SNPDict = {}
    with open(SNPFile) as f1:
        next(f1)
        for variations in f1:
            variationCount += 1
            try:
                chromosome, start, stop, name, ref, altCount, alts, shiftbases, freqSourceCount, minorAlleleFreq, MajorAllele, minorAllele, maxFuncImpact, Class, ucscNotes, dataOffset, dataLen = variations.rstrip().split("\t") 

            except ValueError:
                print(variations)
                break            
            chromosome  = chromosome
            start       = start
            stop        = stop

            if Class == "snv":
                SNPlocation = start #check dit even met MOkry of dit klopt?
                SNPcount += 1
                if chromosome not in SNPDict:
                    SNPDict[chromosome] = []

                SNPDict[chromosome].append(SNPlocation)
           
            if variationCount % 1000000 == 0:
                print("processed {} variations".format(variationCount))
    print("snpDict generated")
    return SNPDict

def generateVariantsdict(variantFile):
    tempDict = {}
    with open(variantFile, "r") as f1:
        for lines in f1:
            chrom, start, stop = lines.rstrip().split("\t")

            if chrom not in tempDict:
                tempDict[chrom] = {}


            tempDict[chrom][start] = []
            #tempDict[chrom].append(start)
    print("variantDict generated")
    return tempDict

def generateTranslationDict(translationFile):
    translationDictionary = {}
    with open(translationFile, "r") as t1:
        next(t1)
        for lines in t1:
            metaData = lines.rstrip().split(",")
            cellInfo = metaData[0].strip('"').split(".")
            variantFileFormat = "sorted_run{}_cell{}_variants".format(cellInfo[1][1:],cellInfo[2] )
            metaData = "{}.{}.{}".format(metaData[-1].replace(" ", "_"),cellInfo[2], cellInfo[1][1:])
            translationDictionary[variantFileFormat] = metaData
    
    return translationDictionary

def generateGeneTranslationDict(refseqFile):
    geneDict = {}
    with open(refseqFile, "r") as f1:
        for lines in f1:
            data = lines.rstrip().split("\t")
            gene = data[-4]
            chrom =data[2]
            start= int(data[4])
            stop = int(data[5])

            if chrom not in geneDict:
                geneDict[chrom] = {}

            geneDict[chrom][start] = {"stop" : stop, "gene" : gene}
    return geneDict

def checkSamples2(variantFolder, variantDict, transDict, snpDict):
    allFiles = os.listdir(variantFolder)
    for variantFile in allFiles:
        #sampleName = "_".join(variantFile.rstrip().split("_")[1:4])
        sampleName = variantFile.split(".")[0]
        with open(variantFolder + "/" + variantFile) as f1:
            next(f1) 
            print(sampleName)           
            for lines in f1:
                chrom, loc, ref, som, A, T, C, G =  lines.rstrip().split("\t")
                #print(loc)
                if chrom in variantDict:
                    if loc in variantDict[chrom]:
                        if sampleName in transDict:
                            variantDict[chrom][loc].append(transDict[sampleName])
                        else:
                            tempSplit = sampleName.split("_")
                            #print(tempSplit)
                            cellRun = "cell{}.{}".format(tempSplit[2][4:], tempSplit[1][3:])
                            variantDict[chrom][loc].append(cellRun)
                        #print("wel")
                    else:
                        pass
                else:
                    pass

        
    with open("cellsPerVariant.txt", "w") as w2:
        for chromosomes in variantDict:
            for variants in variantDict[chromosomes]:
                #if len(variantDict[chromosomes][variants]) == 0:
                #    w2.write("{}_{}\n".format(chromosomes, variants))
                if len(variantDict[chromosomes][variants]) > 1:
                    try:
                        if variants not in snpDict[chromosomes]:
                            w2.write("{}_{}\t{}\t{}\n".format(chromosomes, variants, len(variantDict[chromosomes][variants]) ,"\t".join(variantDict[chromosomes][variants])))
                        else:
                            print("chromosome {} position {} is a common SNP".format(chromosomes, variants))
                    except KeyError:
                            w2.write("{}_{}\t{}\t{}\n".format(chromosomes, variants, len(variantDict[chromosomes][variants]) ,"\t".join(variantDict[chromosomes][variants])))                                        
                else:
                    pass

#example: python3 2_6_variantsPerCell.py db_common_SNPs153.txt AE4450_voor_Mark_somatic_mutations.txt Patient1Variants Patient1_variantsInExons.bed 

#input
commonSNPs          = sys.argv[1] #file with common known SNPS in the format: #chrom	chromStart	chromEnd	name	ref	altCount	alts	shiftBases	freqSourceCount	minorAlleleFreq	majorAllele	minorAllele	maxFuncImpact	class	ucscNotes	_dataOffset	_dataLen
patientCellTypes    = sys.argv[2] #metainformation on cell type
variantsFolder      = sys.argv[3] #folder with all variants
variantsInsideExons = sys.argv[4] #output of python3 2_5_inExons.py

snpDict = getSNPs(commonSNPs) #"db_common_SNPs153.txt"
translationDict = generateTranslationDict(patientCellTypes) #"Patient2_metadata.csv" format: "","orig.ident","nCount_RNA","nFeature_RNA","Patient","Plate","Phenotype","Sex","CD68.score","CD3.score","CD34.score","alpha.SMA.score","SR.score","Glyc.c.score","Calcification","Batch","File","C.H","Type","ID","percent.mito","percent.RPL","percent.RPS","nCount_SCT","nFeature_SCT","SCT_snn_res.0.8","seurat_clusters","new.ident"
variantDict     = generateVariantsdict(variantsInsideExons) #"variantsInExons.bed"  format: #chrom	chromStart	chromEnd	name	ref	altCount	alts	shiftBases	freqSourceCount	minorAlleleFreq	majorAllele	minorAllele	maxFuncImpact	class	ucscNotes	_dataOffset	_dataLen
checkSamples2(variantsFolder, variantDict, translationDict, snpDict)
