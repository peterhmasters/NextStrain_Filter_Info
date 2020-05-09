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

import pandas as pd


class TrackFilter:
    """This class tracks which sequences are included or excluded by the filter and the reason why

    :param seqs: A list of sequences to be filtered on
    :type seqs: list

    :Example:

    >>> seqs = ['PRVABC59','COL/FLR_00008/2015','Colombia/2016/ZC204Se']
    >>> test = TrackFilter(seqs)
    >>> assert test.seqs == seqs
    >>> assert test.df.shape == (len(seqs), 3)


    """

    def __init__(self, seqs):
        """Constructor method

        :keyword df tracking sequence name, status to include or not, and the reason for either
        """
        self.seqs = seqs
        # Initializing each sequenced as included and passed all tests
        self.df = pd.DataFrame(
            list(zip(seqs, ["inclusion"] * len(seqs), ["passed"] * len(seqs))),
            columns=["Sequence", "Status", "Reason"],
        )

    def update(self, action, reason, new_seqs):
        """This function updates the dataframe tracking sequences based on if they need to be excluded or included

        :param action: If these are the sequences that to not be removed ('remove'), or the sequenced to be added ("add")
        :type action: str
        :param reason: The reason why are these sequences being added or removed
        :type reason: str
        :param new_seqs: The list of sequenced that need to not be removed or to be added
        :type new_seqs: list
        :return: None
        :rtype: None

        :Example:

        >>> seqs = ['PRVABC59','COL/FLR_00008/2015','Colombia/2016/ZC204Se']
        >>> test = TrackFilter(seqs)
        >>> test.update('remove', 'too short', ['PRVABC59','COL/FLR_00008/2015'])
        >>> assert test.seqs == ['PRVABC59','COL/FLR_00008/2015']


        """
        # Updating the list of sequences, the self.seqs update isn't used anywhere but is future-proofing
        # In case this class could be used to track sequences instead of global variables
        self.seqs = (
            new_seqs
            if action == "remove"
            else self.seqs + list(filter(lambda seq: seq not in self.seqs, new_seqs))
        )
        # Updating the dataframe if the sequence needs to be removed or it needs to be added and the reason why
        self.df = pd.DataFrame(
            map(
                lambda sequence: [
                    sequence[0],
                    "exclusion" if action == "remove" else "inclusion",
                    reason
                    # only updating the Reason if the Status is changing
                    if ((action == "remove") == (sequence[1] == "inclusion"))
                    else sequence[2],
                ]
                # only updating the sequence if it's remove & not in the sequence or add & in the sequence
                if ((action == "remove") != (sequence[0] in new_seqs)) else sequence,
                self.df.itertuples(index=False),
            ),
            columns=["Sequence", "Status", "Reason"],
        )

    def summary(self, text_output, csv_output):
        """Conditionally printing out a text summary and a csv containing filter information

        :param text_output: To print out the text report or not
        :type text_output: Bool
        :param csv_output: To output the csv report or not
        :type csv_output: Bool
        :return: None
        :rtype: None
        """
        if csv_output:
            self.df.sort_values(["Status", "Reason"]).to_csv("Filter_Results.csv")
        if text_output:
            print(self.df.groupby(["Status", "Reason"]).count())
