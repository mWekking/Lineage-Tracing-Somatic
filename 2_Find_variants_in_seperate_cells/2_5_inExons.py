import sys

#example: python3 2_5_inExons.py patient1_variantWithCellcounts.txt hg19_exons.bed patient1_variantsinExons.bed

#input
variantsFile = sys.argv[1]
exonFile     = sys.argv[2]
outputFile   = sys.argv[3]

exonDict = {}

with open(exonFile, "r") as f1:
    for lines in f1:
        label, chrom, strand, start, stop = lines.rstrip().split(",")

        if "_" in chrom:
            pass
        else:
            if chrom not in exonDict:
                exonDict[chrom] = []

            exonDict[chrom].append([int(start), int(stop)])


with open(variantsFile, "r") as f2, open(outputFile, "w") as w1:
    for lines in f2:
        try:
            chrom, loc = lines.rstrip().split("\t")[0].split("_")
            loc = int(loc)
        except ValueError:
            pass

        for exon in exonDict[chrom]:
            if exon[0] < loc and exon[1] > loc:
                w1.write("{}\t{}\t{}\n".format(chrom, loc, loc + 1))
                break