rule evaluation1:
    input:
        input_kmerCounting1,
        file_aceIDtaxID,
        file_taxIDLineage
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/aceID-taxID-species.tsv"
    run:
        import os
        import pandas as pd
        import dask.dataframe as dd
        import numpy as np
        import time

        acessionID = os.popen(f"cat {input[0]}|grep '^>'|cut -f1 -d' '|sed 's/>//g'").read().split("\n")
        acessionID = list(map(lambda x : x.split(".")[0], acessionID))

        df_taxID = dd.read_csv(input[1], sep="\t", header=0, usecols=["accession", "taxid"], dtype={"accession":"object", "taxid":"int64"})
        df_taxLineage = dd.read_csv(input[2], sep=",", header=0, assume_missing=True, usecols=["tax_id", "species", "no rank"], dtype={"tax_id":"int64", "species":"object", "no rank":"object"})
        df_taxID = df_taxID[df_taxID["accession"].isin(acessionID)]
        df_out = df_taxID.merge(df_taxLineage, how="left", left_on="taxid", right_on="tax_id")

        df_out.compute().to_csv(output[0], sep="\t",index=False)

rule evaluation2:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allOligo.set",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allOligo.set.position",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/aceID-taxID-species.tsv",
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/evaluation/allOligo.set.coverage"
    run:
        import os
        import sys
        import re
        from collections import defaultdict
        from Bio import SeqIO
        from Bio.Seq import Seq
        import pandas as pd
        import numpy as np

        df_oligo = pd.read_table(input[0], sep="\t", header=0, index_col=0)
        df_oligoPosition = pd.read_table(input[1], sep="\t", header=0, index_col=0)
        df_aceIDSpecies = pd.read_table(input[2], sep="\t", header=0, index_col=0)

        df_oligoPosition = df_oligoPosition.sort_values("START")

        dict_aceID_species = df_aceIDSpecies["species"].to_dict()
        dict_species_count = df_aceIDSpecies["species"].value_counts().to_dict()



        # for start, group_df in df_oligoPosition.groupby("START"):
        #     ls_refs = group_df.RNAME.to_list()
