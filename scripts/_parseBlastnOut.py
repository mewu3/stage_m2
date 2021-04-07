import sys

inputFile = sys.argv[1]
outputFile = sys.argv[2]

inputOpen = open(inputFile, "r")
outputOpen = open(outputFile, "w")

pIDSet = set()
for line in inputOpen:
    line = line.rstrip("\n")
    list = line.split("\t")
    pID = list[0]
    print(pID)
    # pID = int(list[0].split("|")[0].lstrip("p"))
    identity = float(list[3])
    if identity > 50:
        pIDSet.add(pID)

print(len(pIDSet))
# print(sorted(pIDSet))
