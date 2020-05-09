from Bio import SeqIO
from TrackFilter import TrackFilter
import pandas as pd
from functools import reduce

seqs = list(SeqIO.to_dict(SeqIO.parse("sequences.fasta", "fasta")).keys())
test = TrackFilter(seqs)
assert test.seqs == seqs
assert test.df.shape == (len(seqs), 3)
test.update("remove", "too short", ["PAN/CDC_259359_V1_V3/2015"])
print(test.seqs)
print(test.df)
test.update("add", "changed my mind", ["PAN/CDC_259359_V1_V3/2015"])
print(test.seqs)
print(test.df)

