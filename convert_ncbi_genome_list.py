#!/usr/bin/env python3

import os
import pandas as pd
import sys

assert len(sys.argv) == 3, "Please specify input and output CSV files"

fp_in = sys.argv[1]
assert os.path.exists(fp_in)

fp_out = sys.argv[2]
assert os.path.exists(fp_out) is False

df = pd.read_csv(fp_in, sep=",")
assert "#Organism Name" in df.columns.values
assert "GenBank FTP" in df.columns.values

df = df.reindex(columns=["#Organism Name", "GenBank FTP"])

df["suffix"] = df["GenBank FTP"].apply(lambda s: s.split("/")[-1])
assert df["suffix"].unique().shape[0] == df.shape[0]


def reformat_name(s, chars=[",", ".", "=", "-", "(", ")", "'", "^", '"', "[", "]"]):
    s = s.strip(" ")
    for c in chars:
        s = s.replace(c, "")

    while "  " in s:
        s = s.replace("  ", " ")
    s = s.replace(" ", "_")
    return s


df["#Organism Name"] = df["#Organism Name"].apply(reformat_name)

# Make sure that all organism names are unique
if df["#Organism Name"].unique().shape[0] < df.shape[0]:
    duplicate_names = df["#Organism Name"].value_counts().index.values[
        df["#Organism Name"].value_counts() > 1
    ]
    if len(duplicate_names) > 0:
        for n in duplicate_names:
            ix = df["#Organism Name"] == n
            df.loc[ix, "#Organism Name"] = df.loc[ix].apply(
                lambda r: "{}_{}".format(r["#Organism Name"], r["suffix"]),
                axis=1
            )

df["fasta"] = df.apply(
    lambda r: "{}/{}_genomic.fna.gz".format(r["GenBank FTP"], r["suffix"]),
    axis=1
)
df["gff"] = df.apply(
    lambda r: "{}/{}_genomic.gff.gz".format(r["GenBank FTP"], r["suffix"]),
    axis=1
)

df.reindex(columns=["#Organism Name", "fasta", "gff"]).to_csv(fp_out, sep=",", header=None, index=None)
