import pandas as pd

with open("./data/sources/HPO/genes_to_phenotype.txt", "r") as infile:
    thing = pd.read_csv(
        "./data/sources/HPO/genes_to_phenotype.txt", sep="\t", header=0, dtype=str)

    fff = {symbol: set(group["hpo_id"].values)
           for symbol, group in thing[
        ["gene_symbol", "hpo_id"]
    ].groupby("gene_symbol")}

    print(fff)
