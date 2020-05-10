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


class FilterInfo:
    """This class tracks which sequences are included or excluded by the filter and the reason why

    :param seqs: A list of sequences to be filtered on
    :type seqs: list

    :Example:

    >>> seqs = ['PRVABC59','COL/FLR_00008/2015','Colombia/2016/ZC204Se']
    >>> filter_info = FilterInfo(seqs)
    >>> assert filter_info.seqs == seqs
    >>> assert filter_info.df.shape == (len(seqs), 3)


    """

    def __init__(self, seqs):
        """Constructor method

        :keyword df tracking sequence name, status to include or not, and the reason for either
        """
        self.seqs = seqs
        # Initializing each sequenced as included and passed all tests
        self.df = pd.DataFrame(
            list(
                zip(seqs, ["inclusion"] * len(seqs), ["passed all filters"] * len(seqs))
            ),
            columns=["Sequence", "Status", "Reason"],
        )

    def update(self, action, reason, new_seqs):
        """This function updates the dataframe tracking sequences based on if they need to be excluded or included

        :param action: If these are the sequences that to not be removed ('filter'), or the sequenced to be added ("add")
        :type action: str
        :param reason: The reason why are these sequences being filtered or added
        :type reason: str
        :param new_seqs: The list of sequenced that need to not be filtered or to be added
        :type new_seqs: list, set, str
        :return: None
        :rtype: None

        :Example:

        >>> seqs = ['PRVABC59','COL/FLR_00008/2015','Colombia/2016/ZC204Se']
        >>> filter_info = FilterInfo(seqs)
        >>> filter_info.update('filter', 'too short', ['PRVABC59','COL/FLR_00008/2015'])
        >>> assert filter_info.seqs == ['PRVABC59','COL/FLR_00008/2015']


        """
        # Checking if the input is a set, if so converting it to a list
        if type(new_seqs) in (set, str):
            new_seqs = list(new_seqs) if type(new_seqs) is set else [new_seqs]
        elif type(new_seqs) is list:
            pass
        else:
            raise TypeError("sequence must be a list, set or str")
        # Updating the list of sequences, the self.seqs update isn't used anywhere but is future-proofing
        # In case this class could be used to track sequences instead of global variables
        self.seqs = (
            new_seqs
            if action == "filter"
            else self.seqs + list(filter(lambda seq: seq not in self.seqs, new_seqs))
        )
        # Updating the dataframe if the sequence needs to be removed or it needs to be added and the reason why
        self.df = pd.DataFrame(
            map(
                lambda sequence: [
                    sequence[0],
                    "exclusion" if action == "filter" else "inclusion",
                    reason
                    # only updating the Reason if the Status is changing
                    if ((action == "filter") == (sequence[1] == "inclusion"))
                    else sequence[2],
                ]
                # only updating the sequence if it's remove & not in the sequence or add & in the sequence
                if ((action == "filter") != (sequence[0] in new_seqs)) else sequence,
                self.df.itertuples(index=False),
            ),
            columns=["Sequence", "Status", "Reason"],
        )

    def summary(self, filter_reasons):
        """Printing out a text summary and conditionally exporting a csv containing filter reasoning

        :param filter_reasons: The filter reasoning destination or None if not wanted
        :type text_output: str
        :return: None
        :rtype: None
        """
        print(self.df.groupby(["Status", "Reason"]).count())
        if filter_reasons:
            self.df.sort_values(["Status", "Reason"]).to_csv(filter_reasons)
