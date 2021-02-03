import sys

cellsPerVariantFile = sys.argv[1]


with open(cellsPerVariantFile, "r") as r1, open("mpileupCoordinates.txt", "w") as w1:
    for lines in r1:
        try:
            chrom, loc = lines.split("\t")[0].split("_")
            w1.write("{} {}\n".format(chrom, loc))
        except ValueError:
            print(lines.split("\t"))