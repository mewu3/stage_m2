input = "/datater/wu/data/enterovirus/primer3/forward_2540-3044.calculated"

import pandas as pd
import os

df = pd.read_table(input, sep="\t", header=0)
mean = float(df["Tm"].mean())
std = float(df["Tm"].std())

TmSeuilPlus = mean + std
TmSeuilLess = mean - std

print(df)
df_filtered = df[(df["Tm"] >= TmSeuilLess) & (df["Tm"] <= TmSeuilPlus) & (df["CG%"] >= 40) & (df["CG%"] <= 60) & (df["hairpin-dG"] > -9000) & (df["homodier-dG"] > -9000)]
print(df_filtered)
