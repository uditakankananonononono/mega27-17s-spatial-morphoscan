"""Appendix: dataset-ledger composition by source with example accessions."""
import pandas as pd
df = pd.read_csv("results/accession_ledger.csv")
NL = chr(92)*2 + chr(10)
with open("paper/ledger_table.tex", "w") as f:
    f.write(chr(92)+"begin{longtable}{p{3.4cm}ccp{7.6cm}}"+chr(92)+"toprule"+chr(10)+"Source & Records & Level & Example accessions "+NL+chr(92)+"midrule"+chr(10)+chr(92)+"endfirsthead"+chr(10)+chr(92)+"toprule Source & Records & Level & Examples "+NL+chr(92)+"midrule "+chr(92)+"endhead"+chr(10))
    for src, g in df.groupby("source"):
        lvl = g["level_of_use"].iloc[0]
        ex = "; ".join(g["accession"].head(3)).replace("_", chr(92)+"_")
        s = src.replace("_", chr(92)+"_")
        f.write(f"{s} & {len(g)} & {lvl} & {ex} "+NL)
    f.write(chr(92)+"midrule Total & "+str(len(df))+" & & "+NL)
    f.write(chr(92)+"bottomrule"+chr(10)+chr(92)+"end{longtable}"+chr(10))
print("written")
