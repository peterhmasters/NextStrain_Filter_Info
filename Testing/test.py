# from Bio import SeqIO
from Testing.TrackFilter import FilterInfo
from tempfile import TemporaryDirectory
import os


# seqs = list(SeqIO.to_dict(SeqIO.parse("sequences.fasta", "fasta")).keys())
seqs = ["PAN/CDC_259359_V1_V3/2015", "COL/FLR_00024/2015", "PRVABC59"]
test = FilterInfo(seqs)
assert test.seqs == seqs
assert test.df.shape == (len(seqs), 3)
test.update("filter", "too short", ["PAN/CDC_259359_V1_V3/2015", "COL/FLR_00024/2015"])
assert test.seqs == ["PAN/CDC_259359_V1_V3/2015", "COL/FLR_00024/2015"]
assert test.df.iloc[0:3, 1:3].values.tolist() == [
    ["inclusion", "passed all filters"],
    ["inclusion", "passed all filters"],
    ["exclusion", "too short"],
]
test.update("filter", "too long", {"PAN/CDC_259359_V1_V3/2015"})
assert test.df.iloc[0:3, 1:3].values.tolist() == [
    ["inclusion", "passed all filters"],
    ["exclusion", "too long"],
    ["exclusion", "too short"],
]
test.update(
    "add", "changed my mind", ["PAN/CDC_259359_V1_V3/2015", "COL/FLR_00024/2015"]
)
assert test.df.iloc[0:3, 1:3].values.tolist() == [
    ["inclusion", "passed all filters"],
    ["inclusion", "changed my mind"],
    ["exclusion", "too short"],
]
test.update("filter", "changed my mind again", "PAN/CDC_259359_V1_V3/2015")
assert test.df.iloc[0:3, 1:3].values.tolist() == [
    ["inclusion", "passed all filters"],
    ["exclusion", "changed my mind again"],
    ["exclusion", "too short"],
]

with TemporaryDirectory() as tmp:
    test.summary(None)
    assert len(os.listdir(tmp)) == 0
    test.summary(os.path.join(tmp, "test.csv"))
    assert len(os.listdir(tmp)) == 1
