# I'm working on the issue: https://github.com/nextstrain/augur/issues/424
# I'm building a class locally to be able to update: https://github.com/nextstrain/augur/blob/master/augur/filter.py

# Initial scratch work I did to understand the dataset
# for record in SeqIO.parse("sequences.fasta", "fasta"):
#     print(record.id)
# # sequences and their metadata
# seqs = SeqIO.to_dict(SeqIO.parse("sequences.fasta", 'fasta'))
# print(seqs)
# # sequences names kept after filters
# seq_keep = list(seqs.keys())
# print(seq_keep)
# # sequences names that we start with
# all_seq = seq_keep.copy()
# print(all_seq)
#
# for seq_name in seq_keep: print(seq_name)

# More scratchwork to set up my plan for this asignment
# End goal is to get sequence: name, status: inclusion /exclusion,  reason : passed all/failed x/added back y
# Plan is to create dataframe with every name, status (inclusion), reason (passed all)
# Then after each filter update the dataframe accordingly
# Then print the dataframe or return the dataframe
# Then print summary stats

# Class of track metadata
# Initiate based on all_seq
# init: List in, df out
# Class method of update
# update: list in, update df
# Class method of print summary stats/df
# output: arg in, print / to_csv / nothing conditional

from Bio import SeqIO
import pandas as pd

class TrackFilter:
    def __init__(self, seqs):
        self.seqs = seqs
        self.df = pd.DataFrame(
            list(zip(seqs, ["inclusion"] * len(seqs), ["passed"] * len(seqs))),
            columns=["Sequence", "Status", "Reason"],
        )


seqs = list(SeqIO.to_dict(SeqIO.parse("sequences.fasta", "fasta")).keys())
test = TrackFilter(seqs)
print(test.seqs)
print(test.df)
