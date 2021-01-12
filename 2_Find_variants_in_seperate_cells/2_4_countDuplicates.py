import sys

#example: python3 2_4_countDuplicates.py run1_cell123_master.txt 2 patient1_variantWithCellcounts.txt

#input
MasterFile      =     sys.argv[1] 
numberOfCells   = int(sys.argv[2])
outputFile      =     sys.argv[3]


scores = []
doneLines = 0


with open(MasterFile, "r") as f1:
    header = f1.readline().rstrip().split("\t")[1:]
    scores = f1.readline().rstrip().split("\t")[1:]
    for lines in f1: 
            print(doneLines)
            addList = lines.rstrip().split("\t")[1:]
            scores = [int(x) + int(y) for x, y in zip(scores, addList)]
            doneLines += 1

with open(outputFile, "w") as f1:
    for i, y in enumerate(scores):
        if y >= numberOfCells:
            f1.write("{}\t-\t{}\n".format(header[i], y))