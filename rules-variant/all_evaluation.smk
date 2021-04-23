rule evaluation:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allOligo.set"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allOligo.set.evaluation"
    run:
        from Bio import SeqIO
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from collections import defaultdict
        import pandas as pd
        import re
