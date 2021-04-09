import sys
import pandas as pd
from collections import defaultdict

input1 = sys.argv[1]
output1 = sys.argv[2]
table = False

if table:
    df = pd.read_table(input1, sep="\t", header=0)
    df = df.sort_values(["position", "kmerCount"], ascending=False)
    df_first = df.groupby("position").first().reset_index()
    df_first=df_first[["position", "kmerCount"]]
    df_first.to_csv(output1, sep="\t", index=False)
else:
    dict = defaultdict(list)
    input1Open = open(input1, "r")
    output1Open = open(output1, "w")

    for line in input1Open:
        line = line.rstrip("\n")
        if line.startswith(">"):
            ls = line.split("|")
            position = ls[1].split()[1]
            count = int(ls[2].split()[1])
            # print(ls)
            # print(position, count)
            dict[position].append(count)

    output1Open.write(f"position\tkmerCount\n")
    for key in dict:
        list = sorted(dict[key])
        output1Open.write(f"{key}\t{list[-1]}\n")

    input1Open.close()
    output1Open.close()
