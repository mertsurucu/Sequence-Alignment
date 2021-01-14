import numpy as np

"""
Align given strings using the Smith-Waterman algorithm.
NB: score matrix and the substitution matrix are different matrices!
"""

def pointers(di, ho, ve):
    pointer = max(di, ho, ve, 0)
    pointer_arr = []
    if di == pointer:
        pointer_arr.append('D')
    if ho == pointer:
        pointer_arr.append('H')
    if ve == pointer:
        pointer_arr.append('V')
    return pointer_arr


class LocalAlignment:
    def __init__(self, string1, string2, gap_penalty, matrix):
        """
        :param string1: first string to be aligned, string
        :param string2: second string to be aligned, string
        :param gap_penalty: gap penalty, integer
        :param matrix: substitution matrix containing scores for amino acid
                       matches and mismatches, dict

        Attention! string1 is used to index columns, string2 is used to index rows
        """
        self.string1 = string1
        self.string2 = string2
        self.gap_penalty = gap_penalty
        self.substitution_matrix = matrix
        self.score_matrix = np.zeros((len(string2) + 1, len(string1) + 1), dtype=np.int)
        self.point_matrix = [[[] for _ in range(len(self.string1) + 1)] for _ in range(len(self.string2) + 1)]
        self.align()

    def align(self):
        """
        Align given strings using the Smith-Waterman algorithm.
        NB: score matrix and the substitution matrix are different matrices!
        """
        self.score_matrix[0, 0:] = [i * self.gap_penalty if self.gap_penalty >= 0 else 0 for i in
                                    range(len(self.string1) + 1)]
        self.score_matrix[0:, 0] = [i * self.gap_penalty if self.gap_penalty >= 0 else 0 for i in
                                    range(len(self.string2) + 1)]
        for f in range(1, len(self.string1) + 1):
            for s in range(1, len(self.string2) + 1):
                di, ve, ho = self.score_matrix[s - 1][f - 1] + self.substitution_matrix[self.string1[f - 1]][
                    self.string2[s - 1]], \
                             self.score_matrix[s - 1][f] + self.gap_penalty, \
                             self.score_matrix[s][f - 1] + self.gap_penalty
                self.score_matrix[s, f] = max(di, ho, ve, 0)
                point_arr = pointers(di, ho, ve)
                self.point_matrix[s][f] = point_arr

    def has_alignment(self):
        """
        :return: True if a local alignment has been found, False otherwise
        """

        return False if np.sum(self.score_matrix) == 0 else True

    def _traceback_rec(self, t, r, str1, str2, x, y, s1='', s2=''):
        if self.score_matrix[y][x] > 0:
            arr = t[y][x]

            if 'D' in arr:
                self._traceback_rec(t, r, str1, str2,
                                    x - 1, y - 1, str1[x - 1] + s1, str2[y - 1] + s2)
            if 'H' in arr:
                self._traceback_rec(t, r, str1, str2,
                                    x - 1, y, str1[x - 1] + s1, '-' + s2)
            if 'V' in arr:
                self._traceback_rec(t, r, str1, str2,
                                    x, y - 1, '-' + s1, str2[y - 1] + s2)
        else:
            r.append((s1, s2))

    def get_alignment(self):
        """
        :return: alignment represented as a tuple of aligned strings
        """
        if self.has_alignment():
            idxs = np.where(self.score_matrix == np.amax(self.score_matrix))
            len1 = int(idxs[1])
            len2 = int(idxs[0])
            results = []
            self._traceback_rec(self.point_matrix, results, self.string1[0:len1], self.string2[0:len2], len1, len2)
            return tuple(results[0])
        else:
            return tuple(["", ""])

    def is_residue_aligned(self, string_number, residue_index):
        """
        :param string_number: number of the string (1 for string1, 2 for string2) to check
        :param residue_index: index of the residue to check
        :return: True if the residue with a given index in a given string has been alined
                 False otherwise
        """
        alignments = self.get_alignment()
        str1, str2 = alignments[0], alignments[1]
        if string_number == 1:
            _str = str1
            string = self.string1
        elif string_number == 2:
            _str = str2
            string = self.string2

        _str = _str.replace('-', '')
        start_idx = string.find(_str)
        if start_idx <= residue_index < start_idx + len(_str):
            return True
        else:
            return False
