import numpy as np

"""
Align given strings using the Needleman-Wunsch algorithm,
store the alignments and the score matrix used to compute those alignments.
NB: score matrix and the substitution matrix are different matrices!
"""

def pointers(di, ho, ve):
    pointer = max(di, ho, ve)
    pointer_arr = []
    if di == pointer:
        pointer_arr.append('D')
    if ho == pointer:
        pointer_arr.append('H')
    if ve == pointer:
        pointer_arr.append('V')
    return pointer_arr


class GlobalAlignment:
    def __init__(self, string1, string2, gap_penalty, matrix):
        """
        :param string1: first string to be aligned, string
        :param string2: second string to be aligned, string
        :param gap_penalty: gap penalty, integer
        :param matrix: substitution matrix containing scores for amino acid
                       matches and mismatches, dict

        Attention! string1 is used to index columns verical, string2 is used to index rows horizontal
        """
        self.string1 = string1
        self.string2 = string2
        self.gap_penalty = gap_penalty
        self.substitution_matrix = matrix
        self.score_matrix = np.zeros((len(self.string2) + 1, len(self.string1) + 1), dtype=np.int)
        self.point_matrix = [[[] for _ in range(len(self.string1) + 1)] for _ in range(len(self.string2) + 1)]
        self.align()

    def align(self):
        """
        Align given strings using the Needleman-Wunsch algorithm,
        store the alignments and the score matrix used to compute those alignments.
        NB: score matrix and the substitution matrix are different matrices!
        """

        self.score_matrix[0, 0:] = [i * self.gap_penalty for i in range(len(self.string1) + 1)]
        self.score_matrix[0:, 0] = [i * self.gap_penalty for i in range(len(self.string2)+1)]
        for f in range(1, len(self.string1) + 1):
            for s in range(1, len(self.string2) + 1):
                di, ve, ho = self.score_matrix[s - 1][f - 1] + self.substitution_matrix[self.string1[f - 1]][self.string2[s - 1]], \
                             self.score_matrix[s - 1][f] + self.gap_penalty, \
                             self.score_matrix[s][f - 1] + self.gap_penalty
                self.score_matrix[s, f] = max(di, ho, ve)
                point_arr = pointers(di, ho, ve)
                self.point_matrix[s][f] = point_arr

    def get_best_score(self):
        """
        :return: the highest score for the aligned strings, int

        """

        return self.score_matrix[-1][-1]

    def get_number_of_alignments(self):
        """
        :return: number of found alignments with the best score
        """

        return len(self.get_alignments())

    def _traceback_rec(self, t, r, str1, str2, x, y, s1='', s2=''):
        if x > 0 or y > 0:
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

    def get_alignments(self):
        """
        :return: list of alignments, where each alignment is represented
                 as a tuple of aligned strings
        """

        len1 = len(self.string1)
        len2 = len(self.string2)
        results = []
        self._traceback_rec(self.point_matrix, results, self.string1, self.string2, len1, len2)

        return results

    def get_score_matrix(self):
        """
        :return: matrix built during the alignment process as a list of lists
        """
        return self.score_matrix.tolist()
