"""Integration tests for aligning module sequences."""
from retromol.alignment import (
    ModuleSequence, 
    MultipleSequenceAlignment,
    parse_primary_sequence
)


def test_pairwise_alignment():
    """Test the pairwise alignment of two module sequences."""
    seq1 = ModuleSequence(
        "seq1",
        parse_primary_sequence([
            "polyketide|B7", 
            "polyketide|B2", 
            "polyketide|A2", 
            "polyketide|D7", 
            "polyketide|B2", 
            "polyketide|B2"
        ])
    )
    seq2 = ModuleSequence(
        "seq2",
        parse_primary_sequence([
            "polyketide|B2", 
            "polyketide|A2", 
            "polyketide|D2", 
            "polyketide|B1", 
            "polyketide|B2",
            "peptide|pubchem|1234567",
        ])
    )
    alignment = seq1.optimal_alignment(seq2, 2, 1)
    assert len(alignment.seq1) == len(alignment.seq2)
    assert len(alignment.seq1) == 7


def test_multiple_sequence_alignment():
    """Test the multiple sequence alignment of three module sequences."""
    seq1 = ModuleSequence(
        "seq1",
        parse_primary_sequence([
            "polyketide|B7", 
            "polyketide|B2", 
            "polyketide|A2", 
            "polyketide|D7", 
            "polyketide|B2", 
            "polyketide|B2"
        ])
    )
    seq2 = ModuleSequence(
        "seq2",
        parse_primary_sequence([
            "polyketide|B2", 
            "polyketide|A2", 
            "polyketide|D2", 
            "polyketide|B1", 
            "polyketide|B2",
            "peptide|pubchem|1234567",
        ])
    )
    seq3 = ModuleSequence(
        "seq3",
        parse_primary_sequence([
            "polyketide|A2", 
            "polyketide|D2",
        ])
    )
    alignment = MultipleSequenceAlignment([seq1, seq2, seq3], 2, 1).get_alignment()

    # Check that all sequences have the same length (i.e., the alignment is correct as gaps are added)
    seq_lens = [len(seq.seq) for seq in alignment]
    assert all([seq_len == seq_lens[0] for seq_len in seq_lens])
