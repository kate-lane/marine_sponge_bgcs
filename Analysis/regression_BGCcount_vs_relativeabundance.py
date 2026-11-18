import pandas as pd
from scipy.stats import linregress

infile = "Genome_Stats_BGC_count.csv" #derived from Supplementary Table 3
outfile = "regression_abundance_vs_BGCcount.csv"
minimum = 4 #calculates slope if group contains at least 4
X = "relative_abundance" #relative abundance is a %
Y = "BGC_count" 

df = pd.read_csv(infile, sep=",")
df["phylum"] = df["taxonomy"].str.extract(r"p__([^;]+)") #parse gtdbtk 
df["class"] = df["taxonomy"].str.extract(r"c__([^;]+)")

groups = [("overall", "all MAGs", df)]
for name, g in df.groupby("phylum"):
    groups.append(("phylum", name, g))
for name, g in df.groupby("class"):
    groups.append(("class", name, g))

rows = []
for level, name, g in groups:
    if len(g) < minimum:
        continue
    fit = linregress(g[X], g[Y])
    rows.append([level, name, len(g), fit.rvalue**2, fit.slope, fit.intercept, fit.pvalue])

stats = pd.DataFrame(rows, columns=["rank", "taxon", "MAG_count", "r2", "slope", "intercept", "pvalue"])
stats.to_csv(outfile, index=False)
