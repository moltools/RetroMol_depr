"""
Alignment module.
"""
import typing as ty 

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage

from retromol_sequencing.fingerprint import CompoundClassMapping

class Matrix:
    def __init__(self, nrows: int, ncols: int) -> None:
        """
        Initialize a matrix of size nrows x ncols.

        :param int nrows: Number of rows in matrix.
        :param int ncols: Number of columns in matrix.
        """
        self._matrix = None
        self._nrows = nrows
        self._ncols = ncols

    def __repr__(self) -> str:
        """
        Return a string representation of the matrix.
        
        :return: String representation of the matrix.
        :rtype: str
        """
        flat_matrix = [str(column) for row in self._matrix for column in row]
        padding = len(max(flat_matrix, key=len))

        display = []
        for row in [list(map(str, row)) for row in self._matrix]:
            row = [(' ' * (padding - len(item)) + item) for item in row]
            display.append(' '.join([str(item) for item in row]))

        return '\n'.join(display)

    def build(self, fill: ty.Union[int, float]) -> None:
        """
        Build the matrix with the given fill value.
        
        :param ty.Union[int, float] fill: Value to fill the matrix with.
        :return: None
        :rtype: None
        """
        self._matrix = [
            [fill for _ in range(self._ncols)]
            for _ in range(self._nrows)
        ]

    def transpose(self) -> ty.List[ty.List[ty.Union[int, float]]]:
        """
        Transpose the matrix.

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

    def add(self, row: int, col: int, value: ty.Union[int, float]) -> ty.List[ty.List[ty.Union[int, float]]]:
        """
        Add a value to the matrix at the given row and column.
        
        :param int row: Row to add value to.
        :param int col: Column to add value to.
        :param ty.Union[int, float] value: Value to add to matrix.
        :return: Matrix with value added.
        :rtype: ty.List[ty.List[ty.Union[int, float]]]
        """
        self._matrix[row][col] = value
        return self._matrix

    def get(self, row: int, col: int) -> ty.Union[int, float]:
        """
        Get the value at the given row and column.
        
        :param int row: Row to get value from.
        :param int col: Column to get value from.
        :return: Value at given row and column.
        :rtype: ty.Union[int, float]
        """
        return self._matrix[row][col]

class AlignmentMatrix(Matrix):
    """
    Alignment matrix.
    """
    def __init__(self, ncols: int, nrows: int) -> None:
        """
        Initialize an alignment matrix of size nrows x ncols.
        
        :param int nrows: Number of rows in matrix.
        :param int ncols: Number of columns in matrix.
        """
        super().__init__(ncols, nrows)

class PairwiseScoreMatrix(Matrix):
    """
    Pairwise score matrix.
    """
    def __init__(self, ncols: int, nrows: int) -> None:
        """
        Initialize a pairwise score matrix of size nrows x ncols.

        :param int nrows: Number of rows in matrix.
        :param int ncols: Number of columns in matrix.
        """
        super().__init__(ncols, nrows)

class ScoringMatrix:
    """
    Scoring matrix.
    """
    def __init__(self, scoring_matrix: ty.Optional[str] = None) -> None:
        """
        Initialize a scoring matrix from a string representation.
        
        :param str scoring_matrix: String representation of scoring matrix.
        """
        if scoring_matrix is not None:
            self._scores = self._parse(scoring_matrix)
        else:
            self._scores = self._parse(CompoundClassMapping.get_scoring_matrix())

    def _parse(self, src: str) -> ty.Dict[str, ty.Dict[str, int]]:
        """
        Parse a scoring matrix from a string representation.
        
        :param str src: String representation of scoring matrix.
        :return: Scoring matrix.
        :rtype: ty.Dict[str, ty.Dict[str, int]]
        """
        scores = {}
        lines = src.strip().split("\n")
        header = lines[0].split(",")

        for i, line in enumerate(lines[1:]):
            motif_a = CompoundClassMapping[header[i]]
            scores[motif_a] = {}

            motifs_b = [CompoundClassMapping[motif] for motif in header]
            pairwise_scores = map(int, line.strip().split(",")[1:]) 
            for j, pairwise_score in zip(motifs_b, pairwise_scores):
                scores[motif_a][j] = pairwise_score

        return scores

    def score(self, motif1: CompoundClassMapping, motif2: CompoundClassMapping) -> float:
        """
        Score two motifs.
        
        :param CompoundClassMapping motif1: First motif.
        :param CompoundClassMapping motif2: Second motif.
        :return: Score of two motifs.
        :rtype: float
        """
        return self._scores[motif1][motif2]

Module = ty.Tuple[CompoundClassMapping, ty.Union[int, None]] # Module with tag.

class ModuleSequence:
    def __init__(self, name: str, module_sequence: ty.Union[ty.List[str], ty.List[Module]]) -> None:
        """
        Initialize a module sequence.
        
        :param str name: Name of module sequence.
        :param ty.Union[ty.List[str], ty.List[Module]] module_sequence: Module sequence.
        """
        self.name = name

        if all([isinstance(x, str) for x in module_sequence]):
            self._seq = self.parse(module_sequence)
        else:
            self._seq = module_sequence

    def tag_idx(self) -> None:
        """
        Tag the modules with their original index in the sequence.
        """
        self._seq = [(module, module_idx) for module_idx, (module, _) in enumerate(self._seq)]

    def clear_tags(self) -> None:
        """
        Clear the tags from the modules.
        """
        self._seq = [(module, None) for module, _ in self._seq]

    def insert(self, idx: int, module: Module) -> None:
        """
        Insert a module at the given index.
        """
        self._seq.insert(idx, module)

    def parse(self, module_sequence: ty.List[str], gap: ty.Optional[str] = '-') -> ty.List[Module]:
        """
        Parse a module sequence from a string representation.
        
        :param ty.List[str] module_sequence: Module sequence as a list of strings.
        :param ty.Optional[str] gap: Gap character.
        :return: Module sequence.
        :rtype: ty.List[Module]
        """
        def get_motif_classification(module: str, gap: ty.Optional[str] = '-'):
            if module == gap:
                return CompoundClassMapping["Undefined"]
            else:
                return CompoundClassMapping[module]

        module_list = []
        for module in module_sequence:
            module_list.append((get_motif_classification(module, gap=gap), None))

        return module_list

    def alignment_matrix(self, other: "ModuleSequence", gap_cost: int, end_gap_cost: int) -> AlignmentMatrix: 
        """
        Create an alignment matrix for the given module sequences.
        
        :param ModuleSequence other: Other module sequence.
        :param int gap_cost: Gap cost.
        :param int end_gap_cost: End gap cost.
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
        scoring_matrix = ScoringMatrix()
        for col in range(1, ncols):
            for row in range(1, nrows):

                # Calculate pairwise score..
                cost = scoring_matrix.score(self._seq[row - 1][0], other._seq[col - 1][0])

                # Calculate penalty score.
                if col == len(other._seq) or row == len(self._seq):
                    penalty = end_gap_cost
                else:
                    penalty = gap_cost

                # Calculate final score and fill in.
                mat.add(row, col, max([
                    mat.get(row - 1, col) - penalty,    # Vertical
                    mat.get(row, col - 1) - penalty,    # Horizontal
                    mat.get(row - 1, col - 1) + cost    # Diagonal
                ]))

        return mat

    def traceback(self, other: "ModuleSequence", row: int, col: int, mat: AlignmentMatrix, gap_cost: int, end_gap_cost: int) -> ty.List[Module]:
        """
        Traceback the alignment matrix.
        
        :param ModuleSequence other: Other module sequence.
        :param int row: Current row.
        :param int col: Current column.
        :param AlignmentMatrix mat: Alignment matrix.
        :param int gap_cost: Gap cost.
        :param int end_gap_cost: End gap cost.
        :return: Aligned module sequence.
        :rtype: ty.List[Module]
        """
        # End stage for one seq when seq is depleted of modules.
        if col == 0:
            return self._seq[:row]
        if row == 0:
            return [(CompoundClassMapping.Undefined, None)] * col

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
            self._seq[row - 1][0] == other._seq[col - 1][0] or
            diagonal > max([horizontal, vertical]) or
            (vertical - current != penalty and horizontal - current != penalty)
        ):
            return self.traceback(other, row - 1, col - 1, mat, gap_cost, end_gap_cost) + [self._seq[row - 1]]

        if horizontal > vertical:
            return self.traceback(other, row, col - 1, mat, gap_cost, end_gap_cost) + [(CompoundClassMapping.Undefined, None)]

        if vertical >= horizontal:
            return self.traceback(other, row - 1, col, mat, gap_cost, end_gap_cost) + [self._seq[row - 1]]

    def optimal_alignment(self, other: "ModuleSequence", gap_cost: int, end_gap_cost: int) -> "PairwiseAlignment":
        """
        Calculate the optimal alignment between two module sequences (Needleman-Wunsch).

        :param ModuleSequence other: Other module sequence.
        :param int gap_cost: Gap cost.
        :param int end_gap_cost: End gap cost.
        :return: Optimal alignment.
        :rtype: PairwiseAlignment
        """
        mat = self.alignment_matrix(other, gap_cost, end_gap_cost)
        aligned_self = self.traceback(other, len(self._seq), len(other._seq), mat, gap_cost, end_gap_cost)
        mat.transpose()
        aligned_other = other.traceback(self, len(other._seq), len(self._seq), mat, gap_cost, end_gap_cost)
        alignment_score = mat.get(-1, -1)
        alignment = PairwiseAlignment(self.name, aligned_self, other.name, aligned_other, alignment_score, gap_cost, end_gap_cost)
        return alignment

class PairwiseAlignment:
    def __init__(self, name1: str, aligned1: ty.List[Module], name2: str, aligned2: ty.List[Module], score: int, gap_cost: int, end_gap_cost: int) -> None:
        """
        Initialize a pairwise alignment.
        
        :param str name1: Name of first module sequence.
        :param ty.List[Module] aligned1: Aligned first module sequence.
        :param str name2: Name of second module sequence.
        :param ty.List[Module] aligned2: Aligned second module sequence.
        :param int score: Alignment score.
        :param int gap_cost: Gap cost.
        :param int end_gap_cost: End gap cost.
        """
        self.name_seq1 = name1
        self.seq1 = aligned1
        self.name_seq2 = name2
        self.seq2 = aligned2
        self.score = score
        self.gap = gap_cost
        self.end_gap = end_gap_cost

    def aligned_sequences(self) -> ty.Tuple["ModuleSequence", "ModuleSequence"]:
        """
        Return the aligned sequences.
        
        :return: Aligned sequences.
        :rtype: ty.Tuple[ModuleSequence, ModuleSequence]
        """
        return (ModuleSequence(self.name_seq1, self.seq1), ModuleSequence(self.name_seq2, self.seq2))

    def percentage_identity(self) -> float:
        """
        Calculate the percentage identity of the alignment.
        
        :return: Percentage identity.
        :rtype: float
        """
        zipped_alignment = list(zip(self.seq1, self.seq2))
        same = sum([(ab == ba) for ((ab, _), (ba, _)) in zipped_alignment])
        identity_score = (same / len(zipped_alignment)) * 100
        return identity_score
    
class MultipleSequenceAlignment:
    """
    Multiple sequence alignment.
    """
    def __init__(self, seqs: ty.List[ModuleSequence], gap_cost: int, gap_end_cost: int) -> None:
        """
        Initialize a multiple sequence alignment.
        
        :param ty.List[ModuleSequence] seqs: List of module sequences to align.
        :param int gap_cost: Gap cost.
        :param int gap_end_cost: End gap cost.
        """
        self._seqs = seqs
        self.gap = gap_cost
        self.end = gap_end_cost
        self.msa = self._align()

    def _align(self) -> None:
        """
        Align the module sequences.
        
        Based on:
        'Progressive Sequence Alignment as a Prerequisite to Correct Phylogenetic Trees' by Feng and Doolittle, 1987
        """
        def all_pairwise_scores(seqs1: ty.List[ModuleSequence], seqs2: ty.List[ModuleSequence]) -> ty.List[CompoundClassMapping]:
            """
            Calculate all pairwise scores between two sets of module sequences.
            
            :param ty.List[ModuleSequence] seqs1: First set of module sequences.
            :param ty.List[ModuleSequence] seqs2: Second set of module sequences.
            :return: Pairwise scores.
            :rtype: ty.List[CompoundClassMapping]
            """
            mat = PairwiseScoreMatrix(len(self._seqs), len(self._seqs))
            mat.build(0.0)
            for idx1, seq1 in enumerate(seqs1):
                for idx2, seq2 in enumerate(seqs2):
                    # Skip second calculation since matrix is mirrored around the diagonal:
                    if idx2 < idx1: continue
                    if idx1 == idx2: score = 0.0
                    else:
                        alignment = seq1.optimal_alignment(seq2, self.gap, self.end)
                        score = 100.0 - alignment.percentage_identity()
                    mat.add(idx1, idx2, score)
                    mat.add(idx2, idx1, score)
            return mat

        # Create a dictionary of all records
        if len(self._seqs) == 0: msa = {0: [ModuleSequence("", "")]} # Return empty alignment.
        elif len(self._seqs) == 1: msa = {0: self._seqs} # Return seq if single.
        else:
            # Identification of most closely related pair.
            mat = all_pairwise_scores(self._seqs, self._seqs)
            guide_tree = linkage(pdist(mat._matrix), method="ward")
            msa = {seq_idx: [seq] for seq_idx, seq in enumerate(self._seqs)}

            # Progressive insertion of neutral elements (can create new gaps, but cannot remove existing gaps).
            for pair_idx, pair in enumerate(guide_tree):

                # Every pair in the guide tree can be a new pair that needs seed alignment or it is a leaf connecting to existing alignment.
                j1, j2 = int(pair[0]), int(pair[1])
                new_idx = pair_idx + len(self._seqs)

                if len(msa[j1]) == 1 and len(msa[j2]) == 1:
                    seed1, seed2 = msa[j1][0], msa[j2][0]
                    alignment = seed1.optimal_alignment(seed2, self.gap, self.end)
                    msa[new_idx] = list(alignment.aligned_sequences())

                elif len(msa[j1]) == 1 or len(msa[j2]) == 1:
                    if len(msa[j1]) == 1: leaf, seqs = msa[j1][0], msa[j2]
                    else: leaf, seqs = msa[j2][0], msa[j1]

                    # Tag already aligned sequences with original location for insertion possible gaps after annealing new seq.
                    front_seq = seqs[0]
                    front_seq.tag_idx()
                    back_seq = seqs[-1]
                    back_seq.tag_idx()

                    front_alignment = front_seq.optimal_alignment(leaf, self.gap, self.end)
                    back_alignment = back_seq.optimal_alignment(leaf, self.gap, self.end)
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
                        if m_tag is None: # New insertion.
                            for seq in other_seqs:
                                seq.insert(m_idx, (CompoundClassMapping.Undefined, None))

                    for seq in other_seqs: seq.clear_tags()
                    anchor.clear_tags()

                    if front_score >= back_score: new_msa = [new, anchor] + other_seqs
                    else: new_msa = other_seqs + [anchor, new]

                    msa[new_idx] = new_msa

                elif len(msa[j1]) != 1 and len(msa[j2]) != 1:
                    # First we need to decide which MSA comes on top. To determine this we need to score the top of msa1 with the
                    # bottom of msa2 and vica versa.
                    msa1, msa2 = msa[j1], msa[j2]

                    # Tag already aligned sequences with original location for insertion possible gaps after annealing new seq.
                    msa1_top, msa1_bottom = msa1[0], msa1[-1]
                    msa2_top, msa2_bottom = msa2[0], msa2[-1]
                    msa1_top.tag_idx()
                    msa1_bottom.tag_idx()
                    msa2_top.tag_idx()
                    msa2_bottom.tag_idx()

                    msa1_top_alignment = msa1_bottom.optimal_alignment(msa2_top, self.gap, self.end)
                    msa2_top_alignment = msa2_bottom.optimal_alignment(msa1_top, self.gap, self.end)
                    msa1_top_score = msa1_top_alignment.percentage_identity()
                    msa2_top_score = msa2_top_alignment.percentage_identity()

                    if msa1_top_score >= msa2_top_score:
                        msa1_align_to, msa2_align_to = msa1_top_alignment.aligned_sequences()
                        msa1_other_seqs = msa1[:-1]
                        msa2_other_seqs = msa2[1:]
                    else:
                        msa2_align_to, msa1_align_to = msa2_top_alignment.aligned_sequences()
                        msa1_other_seqs = msa1[1:]
                        msa2_other_seqs = msa2[:-1]

                    for m_idx, (_, m_tag) in enumerate(msa1_align_to._seq):
                        if m_tag is None: # New insertion.
                            for seq in msa1_other_seqs:
                                seq.insert(m_idx, (CompoundClassMapping.Undefined, None))

                    for m_idx, (_, m_tag) in enumerate(msa2_align_to._seq):
                        if m_tag is None: # New insertion.
                            for seq in msa2_other_seqs:
                                seq.insert(m_idx, (CompoundClassMapping.Undefined, None))

                    for seq in msa1_other_seqs: seq.clear_tags()
                    for seq in msa2_other_seqs: seq.clear_tags()
                    msa1_align_to.clear_tags()
                    msa2_align_to.clear_tags()

                    if msa1_top_score >= msa2_top_score:
                        new_msa = msa1_other_seqs + [msa1_align_to, msa2_align_to] + msa2_other_seqs
                    else:
                        new_msa = msa2_other_seqs + [msa2_align_to, msa1_align_to] + msa1_other_seqs

                    msa[new_idx] = new_msa

                else:
                    raise IndexError('aligned sequences from guide tree failed')

                del msa[j1]
                del msa[j2]

        return msa