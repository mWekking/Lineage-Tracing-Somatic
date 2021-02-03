import sys
import os
from scipy.stats import chisquare
import numpy as np
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()

refSeqFile          = sys.argv[1]    #ncbiREFSEQ.txt
metaDataFile        = sys.argv[2]    #Patient2_metadata.csv
cellsPerVariantFile = sys.argv[3]    #sorted_cellsPerVariant.txt
mainFolderLocation = sys.argv[4]

geneDict = {}
with open(refSeqFile, "r") as f1:
    for lines in f1:
        data = lines.rstrip().split("\t")
        gene = data[-4]
        chrom =data[2]
        start= int(data[4])
        stop = int(data[5])

        if chrom not in geneDict:
            geneDict[chrom] = {}

        geneDict[chrom][start] = {"stop" : stop, "gene" : gene}

countDict = {}
with open(metaDataFile, "r") as f1:
    next(f1)
    for line in f1:
        split = line.strip('"').rstrip().split(",")
        cellInfo = split[-1]

        if cellInfo not in countDict:
            countDict[cellInfo] = 1
        else:
            countDict[cellInfo] += 1


def checkPosInPileup(mpileUpFile, chrom, loc):
    with open(mpileUpFile, "r") as r1:
        next(r1)
        for lines in r1:
            split = lines.split("\t")
            if split[0] != chrom:
                continue
            elif int(split[1]) == loc:
                return True
                break
            elif int(split[1]) > loc:
                return False
                break

def checkPosInAllPileups(mainFolder, chrom, loc):
    directoryLoc = mainFolder
    countDir = {}
    for folders in os.listdir(directoryLoc):
        count = 0 
        #print(folders)
        for files in os.listdir(directoryLoc + "/" + folders):
            
            if checkPosInPileup(directoryLoc + "/" + folders + "/" + files, chrom, loc) == True:
                count  += 1
            else:
                pass
        
        countDir[folders] = count

    return countDir

TcellAbsolute   = 0
for keys in [x for x in countDict.keys() if " t cells" in x.lower()]:
    print(keys)
    TcellAbsolute += countDict[keys]
endoAbsolute   = 0
for keys in [x for x in countDict.keys() if "Endothelial" in x]:
    endoAbsolute += countDict[keys]
macroAbsolute   = 0
for keys in [x for x in countDict.keys() if "Macrophages" in x]:
    macroAbsolute += countDict[keys]
bAbsolute   = 0
for keys in [x for x in countDict.keys() if "B Cells" in x]:
    bAbsolute += countDict[keys]
actaAbsolute   = 0
for keys in [x for x in countDict.keys() if "ACTA2+" in x]:
    actaAbsolute += countDict[keys]
    
allChroms = []
with open(cellsPerVariantFile, "r") as f2, open("counted_cells_perVariant", "w") as w1:
    w1.write("loc\tgene\tcellCount\tTcell(%)\tTRelative\tTcoverage\tBcell(%)\tbRelative\tbCoverage\tmacro(%)\tmRelative\tmCoverage\tEndo(%)\teRelative\teCoverage\tSMcell(%)\tsRelative\tsCoverage\tpValues\tcelltypes\n")    
    for lines in f2:
        TcellCount = 0.0
        endoCount  = 0.0
        macroCount = 0.0
        bCount     = 0.0
        actaCount  = 0.0
        totalCount = 0.0

        data       = lines.rstrip().split("\t")
        loc        = data[0]
        chrom, point = loc.split("_")
        allChroms.append(chrom)
        point = int(point)
        try:
            for genes in geneDict[chrom]:
                if genes < point:
                    if geneDict[chrom][genes]["stop"] > point:
                        gene = geneDict[chrom][genes]["gene"]
                        break
                    else:
                        pass
                else:
                    gene = "unkwn"   
        except KeyError:
            gene = "unkwn"  

        coverageDict = checkPosInAllPileups(mainFolderLocation, chrom, point)      

        cellCount  = data[1]
        cells      = data[2:]
        celltypes  = []

        for celltype in cells:
            if "_t_cells" in celltype.lower():
                TcellCount += 1.0
                print(celltype)

            if "endothelial" in celltype.lower():
                endoCount += 1.0

            if "macrophages" in celltype.lower():
                macroCount += 1.0

            if "_b_cells" in celltype.lower():
                bCount += 1.0

            if "acta2+" in celltype.lower():
                actaCount += 1.0

            celltypes.append(celltype)
            totalCount += 1.0

        tCellsRelative = "{}/{}".format(int(TcellCount), TcellAbsolute)
        endoRelative   = "{}/{}".format(int(endoCount) , endoAbsolute)
        macroRelative  = "{}/{}".format(int(macroCount), macroAbsolute)
        bRelative      = "{}/{}".format(int(bCount)    , bAbsolute)
        actaRelative   = "{}/{}".format(int(actaCount) , actaAbsolute)

        tCoverage      = "{}/{}".format(coverageDict["tcells"], TcellAbsolute)
        endoCoverage   = "{}/{}".format(coverageDict["endothelial"], endoAbsolute)
        macroCoverage  = "{}/{}".format(coverageDict["macrophages"], macroAbsolute)
        bCoverage      = "{}/{}".format(coverageDict["bcells"], bAbsolute)
        actaCoverage   = "{}/{}".format(coverageDict["acta2"], actaAbsolute)

        # chiTest = chisquare(f_obs=[int(TcellCount), int(endoCount), int(macroCount), int(bCount), int(actaCount)],
        #                     f_exp=[coverageDict["tcells"] - int(TcellCount), coverageDict["endothelial"] - int(endoCount), coverageDict["macrophages"] - int(macroCount), coverageDict["bcells"] - int(bCount), coverageDict["acta2"] - int(actaCount)])
        
        stats = importr('stats')
        m = np.array([[int(TcellCount),coverageDict["tcells"]       - int(TcellCount)],
                      [int(endoCount) ,coverageDict["endothelial"]  - int(endoCount)],
                      [int(macroCount),coverageDict["macrophages"]  - int(macroCount)],
                      [int(bCount)    ,coverageDict["bcells"]       - int(bCount)],
                      [int(actaCount) ,coverageDict["acta2"]        - int(actaCount)] 
                      ])
        if not np.any(m) == False and (m >= 0).all():
            fischer = stats.fisher_test(m)
            pValue  = "{:.2e}".format(fischer[0][0])
        else:
            print(loc)
            print(m)
            pValue = "NaN"
        
        w1.write("{}\t{}\t{}\t{:.1f}\t{}\t{}\t{:.1f}\t{}\t{}\t{:.1f}\t{}\t{}\t{:.1f}\t{}\t{}\t{:.1f}\t{}\t{}\t{}\t{}\n".format(
                                               loc,
                                               gene, 
                                               cellCount,
                                               TcellCount / totalCount * 100,
                                               tCellsRelative,
                                               tCoverage,
                                               bCount /  totalCount * 100,
                                               bRelative,
                                               bCoverage,
                                               macroCount /  totalCount * 100,
                                               macroRelative,
                                               macroCoverage,
                                               endoCount /  totalCount * 100,
                                               endoRelative,
                                               endoCoverage,
                                               actaCount /  totalCount * 100,
                                               actaRelative,
                                               actaCoverage,
                                               pValue,
                                               ",".join(celltypes)))


print(set(allChroms))