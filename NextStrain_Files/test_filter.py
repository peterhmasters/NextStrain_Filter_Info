import augur.filter

import pytest

import os

from tempfile import TemporaryDirectory


@pytest.fixture
def mock_priorities_file_valid(mocker):
    mocker.patch(
        "builtins.open", mocker.mock_open(read_data="strain1 5\nstrain2 6\nstrain3 8\n")
    )


@pytest.fixture
def mock_priorities_file_malformed(mocker):
    mocker.patch("builtins.open", mocker.mock_open(read_data="strain1 X\n"))


@pytest.fixture
def mock_run_shell_command(mocker):
    mocker.patch("augur.filter.run_shell_command")


class TestFilter:
    def test_read_vcf_compressed(self):
        seq_keep, all_seq = augur.filter.read_vcf(
            "tests/builds/tb/data/lee_2015.vcf.gz"
        )

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"
        assert seq_keep == all_seq

    def test_read_vcf_uncompressed(self):
        seq_keep, all_seq = augur.filter.read_vcf("tests/builds/tb/data/lee_2015.vcf")

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"
        assert seq_keep == all_seq

    def test_read_priority_scores_valid(self, mock_priorities_file_valid):
        # builtins.open is stubbed, but we need a valid file to satisfy the existence check
        priorities = augur.filter.read_priority_scores(
            "tests/builds/tb/data/lee_2015.vcf"
        )

        assert priorities == {"strain1": 5, "strain2": 6, "strain3": 8}
        assert priorities["strain1"] == 5
        assert priorities["strain42"] == 0, "Default priority is 0 for unlisted sequences"

    def test_read_priority_scores_malformed(self, mock_priorities_file_malformed):
        with pytest.raises(ValueError):
            # builtins.open is stubbed, but we need a valid file to satisfy the existence check
            augur.filter.read_priority_scores("tests/builds/tb/data/lee_2015.vcf")

    def test_read_priority_scores_does_not_exist(self):
        with pytest.raises(FileNotFoundError):
            augur.filter.read_priority_scores("/does/not/exist.txt")

    def test_write_vcf_compressed_input(self, mock_run_shell_command):
        augur.filter.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf.gz", "output_file.vcf.gz", []
        )

        augur.filter.run_shell_command.assert_called_once_with(
            "vcftools --gzvcf tests/builds/tb/data/lee_2015.vcf.gz --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_uncompressed_input(self, mock_run_shell_command):
        augur.filter.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf.gz", []
        )

        augur.filter.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_compressed_output(self, mock_run_shell_command):
        augur.filter.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf.gz", []
        )

        augur.filter.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_uncompressed_output(self, mock_run_shell_command):
        augur.filter.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf", []
        )

        augur.filter.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout  > output_file.vcf",
            raise_errors=True,
        )

    def test_write_vcf_dropped_samples(self, mock_run_shell_command):
        augur.filter.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf", ["x", "y", "z"]
        )

        augur.filter.run_shell_command.assert_called_once_with(
            "vcftools --remove-indv x --remove-indv y --remove-indv z --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout  > output_file.vcf",
            raise_errors=True,
        )

    def test_filter_info(self):
        seqs = ["PAN/CDC_259359_V1_V3/2015", "COL/FLR_00024/2015","PRVABC59"]
        test = augur.filter.FilterInfo(seqs)
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
        with self.assertRaises(TypeError):
            test.update("filter", "changed my mind again", 42)

    def test_filter_info_output(self):
        seqs = ["PAN/CDC_259359_V1_V3/2015", "COL/FLR_00024/2015", "PRVABC59"]
        test = augur.filter.FilterInfo(seqs)
        with TemporaryDirectory() as tmp:
            test.summary(None)
            assert len(os.listdir(tmp)) == 0
            test.summary(os.path.join(tmp, "test.csv"))
            assert len(os.listdir(tmp)) == 1
