from Bio import SeqIO
from TrackFilter import TrackFilter

seqs = list(SeqIO.to_dict(SeqIO.parse("sequences.fasta", "fasta")).keys())
test = TrackFilter(seqs)
assert test.seqs == seqs
assert test.df.shape == (len(seqs), 3)
test.update("remove", "too short", ["PAN/CDC_259359_V1_V3/2015","COL/FLR_00024/2015"])
assert test.seqs == ["PAN/CDC_259359_V1_V3/2015","COL/FLR_00024/2015"]
assert test.df.iloc[0,1] == "inclusion"
assert test.df.iloc[0,2] == "passed"
assert test.df.iloc[1,1] == "inclusion"
assert test.df.iloc[1,2] == "passed"
assert test.df.iloc[2,1] == "exclusion"
assert test.df.iloc[2,2] == "too short"
test.update("remove", "too long", ["PAN/CDC_259359_V1_V3/2015"])
assert test.df.iloc[0,1] == "inclusion"
assert test.df.iloc[0,2] == "passed"
assert test.df.iloc[1,1] == "exclusion"
assert test.df.iloc[1,2] == "too long"
assert test.df.iloc[2,1] == "exclusion"
assert test.df.iloc[2,2] == "too short"
test.update("add", "changed my mind", ["PAN/CDC_259359_V1_V3/2015","COL/FLR_00024/2015"])
assert test.df.iloc[0,1] == "inclusion"
assert test.df.iloc[0,2] == "passed"
assert test.df.iloc[1,1] == "inclusion"
assert test.df.iloc[1,2] == "changed my mind"
assert test.df.iloc[2,1] == "exclusion"
assert test.df.iloc[2,2] == "too short"
test.summary(True, True)
