"""This module contains functions for sequence alignment."""

import typing as ty
from dataclasses import dataclass
from enum import Enum, auto

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage


class Motif:
    """Base class for motifs."""

    ...


class PolyketideType(Enum):
    """Enum for polyketide types."""

    A = auto()
    B = auto()
    C = auto()
    D = auto()
    Undefined = auto()


class PeptideType(Enum):
    """Enum for peptide types."""

    PolarAndCharged = auto()
    SmallHydrophobic = auto()
    SmallNonHydrophobic = auto()
    Tiny = auto()
    Bulky = auto()
    CyclicAliphatic = auto()
    CysteineLike = auto()
    Undefined = auto()

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name


@dataclass
class PolyketideMotif(Motif):
    """Dataclass for polyketide motifs."""

    type: PolyketideType
    decoration_type: int | None

    @classmethod
    def from_dict(cls, data: ty.Dict[str, ty.Any]) -> "PolyketideMotif":
        """Create a PolyketideMotif object from a dictionary.

        :param data: Dictionary with polyketide motif data.
        :type data: ty.Dict[str, ty.Any]
        :return: Polyketide motif object.
        :rtype: PolyketideMotif
        """
        accesory_domains = data["accessory_domains"]

        # Determine pk_type.
        if len(accesory_domains) == 0:
            pk_type = PolyketideType.A

        elif len(accesory_domains) == 1 and "KR" in accesory_domains:
            pk_type = PolyketideType.B

        elif (
            len(accesory_domains) == 2
            and "KR" in accesory_domains
            and "DH" in accesory_domains
        ):
            pk_type = PolyketideType.C

        elif (
            len(accesory_domains) == 3
            and "KR" in accesory_domains
            and "DH" in accesory_domains
            and "ER" in accesory_domains
        ):
            pk_type = PolyketideType.D

        else:
            pk_type = PolyketideType.Undefined

        # Determine decoration_type.
        try:
            decoration_type = int(data["decoration_type"])
        except ValueError:
            decoration_type = None

        return cls(pk_type, decoration_type)

    def __str__(self) -> str:
        return f"{self.type}{self.decoration_type}"

    def __repr__(self) -> str:
        return f"{self.type}{self.decoration_type}"


@dataclass
class PeptideMotif(Motif):
    """Dataclass for peptide motifs."""

    type: PeptideType | None

    @classmethod
    def from_dict(cls, data: ty.Dict[str, ty.Any]) -> "PeptideMotif":
        """Create a PeptideMotif object from a dictionary.

        :param data: Dictionary with peptide motif data.
        :type data: ty.Dict[str, ty.Any]
        :return: Peptide motif object.
        :rtype: PeptideMotif
        """
        classification = data["classification"]

        if classification == "small non-hydrophobic":
            peptide_type = PeptideType.SmallNonHydrophobic

        elif classification == "small hydrophobic":
            peptide_type = PeptideType.SmallHydrophobic

        elif classification == "cyclic aliphatic":
            peptide_type = PeptideType.CyclicAliphatic

        elif classification == "cysteine-like":
            peptide_type = PeptideType.CysteineLike

        elif classification == "polar and charged":
            peptide_type = PeptideType.PolarAndCharged

        elif classification == "tiny":
            peptide_type = PeptideType.Tiny

        elif classification == "bulky":
            peptide_type = PeptideType.Bulky

        else:
            peptide_type = PeptideType.Undefined

        return cls(peptide_type)

    def __str__(self) -> str:
        return f"{self.type}"

    def __repr__(self) -> str:
        return f"{self.type}"


class UndefinedMotif(Motif):
    """Class for undefined motifs."""

    ...


class Gap(Motif):
    """Class for gaps."""

    ...


def score_motif_similarity(motif1: Motif, motif2: Motif) -> int:
    """Score the similarity between two motifs.

    :param motif1: First motif.
    :type motif1: Motif
    :param motif2: Second motif.
    :type motif2: Motif
    :return: Score for similarity.
    :rtype: int
    """
    # Both are gaps.
    if isinstance(motif1, Gap) and isinstance(motif2, Gap):
        return -2

    # One is a gap.
    if isinstance(motif1, Gap) or isinstance(motif2, Gap):
        return 0

    # Both are any polyketide motif.
    if isinstance(motif1, PolyketideMotif) and isinstance(motif2, PolyketideMotif):
        if (
            motif1.type == motif2.type
            and motif1.decoration_type == motif2.decoration_type
        ):
            return 3
        elif (
            motif1.type == motif2.type
            and motif1.decoration_type != motif2.decoration_type
        ):
            return 2
        else:
            return 1

    # Both are any peptide motif.
    if isinstance(motif1, PeptideMotif) and isinstance(motif2, PeptideMotif):
        if motif1.type == motif2.type:
            return 3
        else:
            return 1

    # One is a polyketide motif and the other is a peptide motif.
    if (
        isinstance(motif1, PolyketideMotif)
        and isinstance(motif2, PeptideMotif)
        or isinstance(motif1, PeptideMotif)
        and isinstance(motif2, PolyketideMotif)
    ):
        return -1

    return 0


def parse_primary_sequence(sequence: ty.Dict[str, ty.Any]) -> "ModuleSequence":
    """
    Parse primary sequence from dictionary.

    :param sequence: Dictionary with primary sequence data.
    :type sequence: ty.Dict[str, ty.Any]
    :return: ModuleSequence object.
    :rtype: ModuleSequence
    """
    module_list = []

    for props in sequence:
        if identifier := props.get("identifier"):
            if identifier == "polyketide":
                module_list.append((PolyketideMotif.from_dict(props), None))
            elif identifier == "peptide":
                module_list.append((PeptideMotif.from_dict(props), None))
            else:
                module_list.append((UndefinedMotif(), None))
        else:
            module_list.append((UndefinedMotif(), None))

    return module_list


class Matrix:
    """Base class for matrices.

    :param nrows: Number of rows.
    :type nrows: int
    :param ncols: Number of columns.
    :type ncols: int
    """

    def __init__(self, nrows: int, ncols: int) -> None:
        self._matrix = None
        self._nrows = nrows
        self._ncols = ncols

    def build(self, fill: ty.Union[int, float]) -> None:
        """Build the matrix with a fill value.

        :param fill: Fill value.
        :type fill: ty.Union[int, float]
        :return: None
        :rtype: None
        """
        self._matrix = [[fill for _ in range(self._ncols)] for _ in range(self._nrows)]

    def transpose(self) -> ty.List[ty.List[ty.Union[int, float]]]:
        """Transpose the matrix.

        :return: Transposed matrix.
        :rtype: ty.List[ty.List[ty.Union[int, float]]]
        """
        transposed_matrix = [[] for _ in range(self._ncols)]
        for row in self._matrix:
            for i, column in enumerate(row):
                transposed_matrix[i].append(column)

        self._nrows, self._ncols = self._ncols, self._nrows
        self._matrix = transposed_matrix

        return self._matrix

    def add(
        self, row: int, col: int, value: ty.Union[int, float]
    ) -> ty.List[ty.List[ty.Union[int, float]]]:
        """Add a value to the matrix.

        :param row: Row index.
        :type row: int
        :param col: Column index.
        :type col: int
        :param value: Value to add.
        :type value: ty.Union[int, float]
        :return: Matrix with added value.
        :rtype: ty.List[ty.List[ty.Union[int, float]]]
        """
        self._matrix[row][col] = value
        return self._matrix

    def get(self, row: int, col: int) -> ty.Union[int, float]:
        """Get a value from the matrix.

        :param row: Row index.
        :type row: int
        :param col: Column index.
        :type col: int
        :return: Value from the matrix.
        :rtype: ty.Union[int, float]
        """
        return self._matrix[row][col]


class AlignmentMatrix(Matrix):
    """Matrix for sequence alignment.

    :param nrows: Number of rows.
    :type nrows: int
    :param ncols: Number of columns.
    :type ncols: int
    """

    def __init__(self, ncols: int, nrows: int) -> None:
        super().__init__(ncols, nrows)


class PairwiseScoreMatrix(Matrix):
    """Matrix for pairwise scores.

    :param ncols: Number of columns.
    :type ncols: int
    :param nrows: Number of rows.
    :type nrows: int
    """

    def __init__(self, ncols: int, nrows: int) -> None:
        super().__init__(ncols, nrows)


Module = ty.Tuple[Motif, ty.Union[int, None]]  # Module with tag.


class ModuleSequence:
    """Class for module sequences.

    :param name: Name of the sequence.
    :type name: str
    :param module_sequence: Sequence of modules.
    :type module_sequence: ty.List[Module]
    """

    def __init__(self, name: str, module_sequence: ty.List[Module]) -> None:
        self.name = name
        self._seq = module_sequence

    def tag_idx(self) -> None:
        """Tag the modules with their original index in the sequence.

        :return: None
        :rtype: None
        """
        self._seq = [
            (module, module_idx) for module_idx, (module, _) in enumerate(self._seq)
        ]

    def clear_tags(self) -> None:
        """Clear the tags from the modules.

        :return: None
        :rtype: None
        """
        self._seq = [(module, None) for module, _ in self._seq]

    def insert(self, idx: int, module: Module) -> None:
        """Insert a module at a specific index.

        :param idx: Index to insert the module.
        :type idx: int
        :param module: Module to insert.
        :type module: Module
        :return: None
        :rtype: None
        """
        self._seq.insert(idx, module)

    def alignment_matrix(
        self, other: "ModuleSequence", gap_cost: int, end_gap_cost: int
    ) -> AlignmentMatrix:
        """Create an alignment matrix.

        :param other: Other sequence.
        :type other: ModuleSequence
        :param gap_cost: Gap cost.
        :type gap_cost: int
        :param end_gap_cost: End gap cost.
        :type end_gap_cost: int
        :return: Alignment matrix.
        :rtype: AlignmentMatrix
        """
        # Instantiate zero matrix.
        nrows = len(self._seq) + 1
        ncols = len(other._seq) + 1
        mat = AlignmentMatrix(nrows, ncols)
        mat.build(0)

        # Fill in initial alignment for when there is zero alignment.
        for row in range(1, nrows):
            mat.add(row, 0, (mat.get(row - 1, 0) - end_gap_cost))

        for col in range(1, ncols):
            mat.add(0, col, (mat.get(0, col - 1) - end_gap_cost))

        # Go through matrix consecutively and calculate the scores.
        for col in range(1, ncols):
            for row in range(1, nrows):

                # Calculate pairwise score..
                cost = score_motif_similarity(
                    self._seq[row - 1][0], other._seq[col - 1][0]
                )

                # Calculate penalty score.
                if col == len(other._seq) or row == len(self._seq):
                    penalty = end_gap_cost
                else:
                    penalty = gap_cost

                # Calculate final score and fill in.
                mat.add(
                    row,
                    col,
                    max(
                        [
                            mat.get(row - 1, col) - penalty,  # Vertical
                            mat.get(row, col - 1) - penalty,  # Horizontal
                            mat.get(row - 1, col - 1) + cost,  # Diagonal
                        ]
                    ),
                )

        return mat

    def traceback(
        self,
        other: "ModuleSequence",
        row: int,
        col: int,
        mat: AlignmentMatrix,
        gap_cost: int,
        end_gap_cost: int,
    ) -> ty.List[Module]:
        """Traceback through the alignment matrix.

        :param other: Other sequence.
        :type other: ModuleSequence
        :param row: Row index.
        :type row: int
        :param col: Column index.
        :type col: int
        :param mat: Alignment matrix.
        :type mat: AlignmentMatrix
        :param gap_cost: Gap cost.
        :type gap_cost: int
        :param end_gap_cost: End gap cost.
        :type end_gap_cost: int
        :return: Aligned sequence.
        :rtype: ty.List[Module]
        """
        # End stage for one seq when seq is depleted of modules.
        if col == 0:
            return self._seq[:row]
        if row == 0:
            return [(Gap, None)] * col

        # Gap penalty depends on location in matrix.
        if row == len(self._seq) or col == len(other._seq):
            penalty = end_gap_cost
        else:
            penalty = gap_cost

        # Define current location and directions.
        current = mat.get(row, col)
        horizontal = mat.get(row, col - 1)
        diagonal = mat.get(row - 1, col - 1)
        vertical = mat.get(row - 1, col)

        # Traceback defined by conditions for every direction.
        if (
            self._seq[row - 1][0] == other._seq[col - 1][0]
            or diagonal > max([horizontal, vertical])
            or (vertical - current != penalty and horizontal - current != penalty)
        ):
            return self.traceback(
                other, row - 1, col - 1, mat, gap_cost, end_gap_cost
            ) + [self._seq[row - 1]]

        if horizontal > vertical:
            return self.traceback(other, row, col - 1, mat, gap_cost, end_gap_cost) + [
                (Gap, None)
            ]

        if vertical >= horizontal:
            return self.traceback(other, row - 1, col, mat, gap_cost, end_gap_cost) + [
                self._seq[row - 1]
            ]

    def optimal_alignment(
        self, other: "ModuleSequence", gap_cost: int, end_gap_cost: int
    ) -> "PairwiseAlignment":
        """Calculate the optimal alignment between two sequences.

        :param other: Other sequence.
        :type other: ModuleSequence
        :param gap_cost: Gap cost.
        :type gap_cost: int
        :param end_gap_cost: End gap cost.
        :type end_gap_cost: int
        :return: Pairwise alignment.
        :rtype: PairwiseAlignment
        """
        mat = self.alignment_matrix(other, gap_cost, end_gap_cost)
        aligned_self = self.traceback(
            other, len(self._seq), len(other._seq), mat, gap_cost, end_gap_cost
        )
        mat.transpose()
        aligned_other = other.traceback(
            self, len(other._seq), len(self._seq), mat, gap_cost, end_gap_cost
        )
        alignment_score = mat.get(-1, -1)
        alignment = PairwiseAlignment(
            self.name,
            aligned_self,
            other.name,
            aligned_other,
            alignment_score,
            gap_cost,
            end_gap_cost,
        )
        return alignment


class PairwiseAlignment:
    """Class for pairwise alignment.

    :param name1: Name of the first sequence.
    :type name1: str
    :param aligned1: Aligned sequence 1.
    :type aligned1: ty.List[Module]
    :param name2: Name of the second sequence.
    :type name2: str
    :param aligned2: Aligned sequence 2.
    :type aligned2: ty.List[Module]
    :param score: Alignment score.
    :type score: int
    :param gap_cost: Gap cost.
    :type gap_cost: int
    :param end_gap_cost: End gap cost.
    :type end_gap_cost: int
    """

    def __init__(
        self,
        name1: str,
        aligned1: ty.List[Module],
        name2: str,
        aligned2: ty.List[Module],
        score: int,
        gap_cost: int,
        end_gap_cost: int,
    ) -> None:
        self.name_seq1 = name1
        self.seq1 = aligned1
        self.name_seq2 = name2
        self.seq2 = aligned2
        self.score = score
        self.gap = gap_cost
        self.end_gap = end_gap_cost

    def aligned_sequences(self) -> ty.Tuple["ModuleSequence", "ModuleSequence"]:
        """Get the aligned sequences.

        :return: Aligned sequences.
        :rtype: ty.Tuple[ModuleSequence, ModuleSequence]
        """
        return (
            ModuleSequence(self.name_seq1, self.seq1),
            ModuleSequence(self.name_seq2, self.seq2),
        )

    def percentage_identity(self) -> float:
        """Calculate the percentage identity of the alignment.

        :return: Percentage identity.
        :rtype: float
        """
        zipped_alignment = list(zip(self.seq1, self.seq2))
        same = sum([(ab == ba) for ((ab, _), (ba, _)) in zipped_alignment])
        return (same / len(zipped_alignment)) * 100


class MultipleSequenceAlignment:
    """Class for multiple sequence alignment.

    :param seqs: List of sequences.
    :type seqs: ty.List[ModuleSequence]
    :param gap_cost: Gap cost.
    :type gap_cost: int
    :param gap_end_cost: End gap cost.
    :type gap_end_cost: int
    """

    def __init__(
        self, seqs: ty.List[ModuleSequence], gap_cost: int, gap_end_cost: int
    ) -> None:
        self._seqs = seqs
        self.gap = gap_cost
        self.end = gap_end_cost
        self.msa = self._align()

    def _align(self) -> None:
        """Align the sequences.

        Based on: 'Progressive Sequence Alignment as a Prerequisite to Correct
        Phylogenetic Trees' by Feng and Doolittle, 1987
        """

        def all_pairwise_scores(
            seqs1: ty.List[ModuleSequence], seqs2: ty.List[ModuleSequence]
        ) -> ty.List[Motif]:
            """
            Calculate all pairwise scores between two sets of sequences.

            :param seqs1: First set of sequences.
            :type seqs1: ty.List[ModuleSequence]
            :param seqs2: Second set of sequences.
            :type seqs2: ty.List[ModuleSequence]
            :return: Pairwise score matrix.
            :rtype: ty.List[Motif]
            """
            mat = PairwiseScoreMatrix(len(self._seqs), len(self._seqs))
            mat.build(0.0)
            for idx1, seq1 in enumerate(seqs1):
                for idx2, seq2 in enumerate(seqs2):
                    # Skip second calculation since matrix is mirrored around
                    # the diagonal:
                    if idx2 < idx1:
                        continue
                    if idx1 == idx2:
                        score = 0.0
                    else:
                        alignment = seq1.optimal_alignment(seq2, self.gap, self.end)
                        score = 100.0 - alignment.percentage_identity()
                    mat.add(idx1, idx2, score)
                    mat.add(idx2, idx1, score)
            return mat

        # Create a dictionary of all records
        if len(self._seqs) == 0:
            msa = {0: [ModuleSequence("", "")]}  # Return empty alignment.
        elif len(self._seqs) == 1:
            msa = {0: self._seqs}  # Return seq if single.
        else:
            # Identification of most closely related pair.
            mat = all_pairwise_scores(self._seqs, self._seqs)
            guide_tree = linkage(pdist(mat._matrix), method="ward")
            msa = {seq_idx: [seq] for seq_idx, seq in enumerate(self._seqs)}

            # Progressive insertion of neutral elements (can create new gaps,
            # but cannot remove existing gaps).
            for pair_idx, pair in enumerate(guide_tree):

                # Every pair in the guide tree can be a new pair that needs
                # seed alignment or it is a leaf connecting
                # to existing alignment.
                j1, j2 = int(pair[0]), int(pair[1])
                new_idx = pair_idx + len(self._seqs)

                if len(msa[j1]) == 1 and len(msa[j2]) == 1:
                    seed1, seed2 = msa[j1][0], msa[j2][0]
                    alignment = seed1.optimal_alignment(seed2, self.gap, self.end)
                    msa[new_idx] = list(alignment.aligned_sequences())

                elif len(msa[j1]) == 1 or len(msa[j2]) == 1:
                    if len(msa[j1]) == 1:
                        leaf, seqs = msa[j1][0], msa[j2]
                    else:
                        leaf, seqs = msa[j2][0], msa[j1]

                    # Tag already aligned sequences with original location for
                    # insertion possible gaps after annealing new seq.
                    front_seq = seqs[0]
                    front_seq.tag_idx()
                    back_seq = seqs[-1]
                    back_seq.tag_idx()

                    front_alignment = front_seq.optimal_alignment(
                        leaf, self.gap, self.end
                    )
                    back_alignment = back_seq.optimal_alignment(
                        leaf, self.gap, self.end
                    )
                    front_score = front_alignment.percentage_identity()
                    back_score = back_alignment.percentage_identity()

                    if front_score >= back_score:
                        align_to = front_alignment
                        other_seqs = seqs[1:]
                    else:
                        align_to = back_alignment
                        other_seqs = seqs[:-1]

                    anchor, new = align_to.aligned_sequences()
                    for m_idx, (_, m_tag) in enumerate(anchor._seq):
                        if m_tag is None:  # New insertion.
                            for seq in other_seqs:
                                seq.insert(m_idx, (Gap, None))

                    for seq in other_seqs:
                        seq.clear_tags()
                    anchor.clear_tags()

                    if front_score >= back_score:
                        new_msa = [new, anchor] + other_seqs
                    else:
                        new_msa = other_seqs + [anchor, new]

                    msa[new_idx] = new_msa

                elif len(msa[j1]) != 1 and len(msa[j2]) != 1:
                    # First we need to decide which MSA comes on top. To
                    # determine this we need to score the top of msa1 with the
                    # bottom of msa2 and vica versa.
                    msa1, msa2 = msa[j1], msa[j2]

                    # Tag already aligned sequences with original location for
                    # insertion possible gaps after annealing
                    # new seq.
                    msa1_top, msa1_bottom = msa1[0], msa1[-1]
                    msa2_top, msa2_bottom = msa2[0], msa2[-1]
                    msa1_top.tag_idx()
                    msa1_bottom.tag_idx()
                    msa2_top.tag_idx()
                    msa2_bottom.tag_idx()

                    msa1_top_alignment = msa1_bottom.optimal_alignment(
                        msa2_top, self.gap, self.end
                    )
                    msa2_top_alignment = msa2_bottom.optimal_alignment(
                        msa1_top, self.gap, self.end
                    )
                    msa1_top_score = msa1_top_alignment.percentage_identity()
                    msa2_top_score = msa2_top_alignment.percentage_identity()

                    if msa1_top_score >= msa2_top_score:
                        msa1_align_to, msa2_align_to = (
                            msa1_top_alignment.aligned_sequences()
                        )
                        msa1_other_seqs = msa1[:-1]
                        msa2_other_seqs = msa2[1:]
                    else:
                        msa2_align_to, msa1_align_to = (
                            msa2_top_alignment.aligned_sequences()
                        )
                        msa1_other_seqs = msa1[1:]
                        msa2_other_seqs = msa2[:-1]

                    for m_idx, (_, m_tag) in enumerate(msa1_align_to._seq):
                        if m_tag is None:  # New insertion.
                            for seq in msa1_other_seqs:
                                seq.insert(m_idx, (Gap, None))

                    for m_idx, (_, m_tag) in enumerate(msa2_align_to._seq):
                        if m_tag is None:  # New insertion.
                            for seq in msa2_other_seqs:
                                seq.insert(m_idx, (Gap, None))

                    for seq in msa1_other_seqs:
                        seq.clear_tags()
                    for seq in msa2_other_seqs:
                        seq.clear_tags()
                    msa1_align_to.clear_tags()
                    msa2_align_to.clear_tags()

                    if msa1_top_score >= msa2_top_score:
                        new_msa = (
                            msa1_other_seqs
                            + [msa1_align_to, msa2_align_to]
                            + msa2_other_seqs
                        )
                    else:
                        new_msa = (
                            msa2_other_seqs
                            + [msa2_align_to, msa1_align_to]
                            + msa1_other_seqs
                        )

                    msa[new_idx] = new_msa

                else:
                    raise IndexError("aligned sequences from guide tree failed")

                del msa[j1]
                del msa[j2]

        return msa

    def get_alignment(self) -> ty.List[ty.List[Module]]:
        """Get the alignment.

        :return: Alignment.
        :rtype: ty.List[ty.List[Module]]
        """
        clusters = list(self.msa.keys())

        if len(clusters) == 1:
            return self.msa[clusters[0]]
        else:
            raise IndexError("Multiple sequence alignment incomplete!")
