"""
Core logic for OOGGA - Optimal Overhang Golden Gate Assembly.

Copyright (C) 2025 Mukundan S <mukundan.kollam@gmail.com>
Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International.
https://creativecommons.org/licenses/by-nc-sa/4.0/

Please cite https://doi.org/10.1101/2025.06.16.659877 if you use this program.
"""


def load_csv_table_as_di(file):
    """Load overhang scoring table from a CSV file.

    Returns a dict mapping overhang sequence -> (relative_efficiency, fidelity).
    """
    with open(file) as fp:
        lines = fp.readlines()

    di = {}
    diagonal_pos = 1
    max_count = 0

    for line in lines[1:]:
        words = line.strip().split(',')
        if max_count < int(words[diagonal_pos]):
            max_count = int(words[diagonal_pos])
        diagonal_pos += 1

    diagonal_pos = 1
    for line in lines[1:]:
        words = line.strip().split(',')
        wc_pair = int(words[diagonal_pos])
        total = sum([int(x) for x in words[1:]])
        relative_efficiency = (wc_pair / max_count) * 100
        fidelity = (wc_pair / total) * 100
        di[words[0]] = (relative_efficiency, fidelity)
        diagonal_pos += 1

    return di


def read_fasta(file):
    """Read a FASTA file and return a dict of {header: sequence}."""
    seqs = {}
    seq = ''
    header = False

    with open(file) as inf:
        for line in inf:
            if line[0] == '>':
                if header:
                    seqs[header] = seq
                header = line.strip()
                seq = ''
            else:
                seq += line.strip()
        seqs[header] = seq.upper()

    return seqs


class Dyna_frag:
    """Dynamic programming fragmenter for Golden Gate assembly.

    Fragments a DNA sequence such that overhangs have optimal
    ligation efficiency and fidelity values.
    """

    def __init__(self, seq, min_len, max_len, overhang_table,
                 n_frag=False, n_trace=5, start_site=0,
                 eff_w=1, fid_w=1, exclusion_list=[],
                 inti_score=1, max_overhang_identity=2,
                 alien_overhangs=[]):
        """
        Parameters
        ----------
        seq : str
            DNA sequence (non-ATGC residues are ignored; overhangs not made there).
        min_len : int
            Minimum fragment length.
        max_len : int
            Maximum fragment length.
        overhang_table : str
            Path to the CSV scoring table (Potapov et al. 2018 supplementary data).
        n_frag : int or False
            Fixed number of fragments. False = auto-determine.
        n_trace : int
            Number of top solutions to report.
        start_site : int
            Index of the first fragment start position.
        eff_w : float
            Exponent weight for efficiency metric.
        fid_w : float
            Exponent weight for fidelity metric.
        exclusion_list : list of int
            DNA indices where overhangs cannot be placed.
        inti_score : float
            Initiation score (default 1).
        max_overhang_identity : int
            Maximum allowed identity between two overhangs.
        alien_overhangs : list of str
            External overhangs to consider when checking identity.
        """
        self.seq = seq
        self.min_len = min_len
        self.max_len = max_len
        self.start_site = start_site
        self.eff_w = eff_w
        self.fid_w = fid_w
        self.n_trace = n_trace
        self.exclusion_list = exclusion_list
        self.inti_score = inti_score
        self.max_overhang_identity = max_overhang_identity
        self._5overhang = ''
        self._3overhang = ''
        self.alien_overhangs = alien_overhangs

        self.score_di = load_csv_table_as_di(overhang_table)
        print('Overhang scores loaded from', overhang_table)

        self.make_matrix()
        print('Matrix populated')

        self.traceback(n_frag)
        print('Trace back completed')

    def make_matrix(self):
        K = self.round_up(len(self.seq) / self.min_len)
        self.mat = [[False for _ in range(len(self.seq))] for _ in range(K + 1)]
        self.split_mat = [[False for _ in range(len(self.seq))] for _ in range(K + 1)]

        self.mat[0][self.start_site] = self.inti_score
        self.split_mat[0][self.start_site] = (self.inti_score, self.inti_score)

        self.trace_di = {}

        for i in range(1, K + 1):
            min_j, max_j = self.get_j_range(i)

            for j in range(len(self.seq)):
                if j < min_j or j > max_j or j in self.exclusion_list:
                    continue

                ovrhg = self.seq[j:j + 4]
                try:
                    eff, fid = self.score_di[ovrhg]
                except KeyError:
                    continue

                sum_scores = sorted(
                    self.__get_score_list(i, j, eff, fid),
                    reverse=True, key=lambda x: x[0]
                )

                if not sum_scores:
                    continue

                score_tot, j_, eff_new, fid_new = sum_scores[0]
                self.mat[i][j] = score_tot
                self.split_mat[i][j] = (eff_new, fid_new)
                self.trace_di[self.make_key(i, j)] = self.make_key(i - 1, j_)

    def __get_score_list(self, i, j, eff, fid):
        score_li = []
        for j_ in range(len(self.seq)):
            frag_len = j - j_
            if (frag_len >= self.min_len and frag_len <= self.max_len
                    and self.mat[i - 1][j_]
                    and self.__overlap_pass(i, j, j_)):
                eff_tally, fid_tally = self.split_mat[i - 1][j_]
                eff_new = eff_tally * (eff / 100)
                fid_new = fid_tally * (fid / 100)
                score = (eff_new ** self.eff_w) * (fid_new ** self.fid_w)
                score_li.append((score, j_, eff_new, fid_new))
        return score_li

    def __overlap_pass(self, i, j, j_):
        js = self.__trace(i - 1, j_)
        current_overhang = self.seq[j:j + 4]
        overhangs = (
            [self.seq[x:x + 4] for x in js]
            + [self.seq[-4:]]
            + self.alien_overhangs
        )
        identities = self.__find_identities(
            overhangs + [self.get_rc(x) for x in overhangs + [current_overhang]],
            current_overhang
        )
        return all(x <= self.max_overhang_identity for x in identities)

    def __find_identities(self, overhangs, current_overhang):
        return [
            sum(1 for x in range(len(current_overhang)) if current_overhang[x] == ovrhg[x])
            for ovrhg in overhangs
        ]

    def get_j_range(self, i):
        scored_prev = [x for x in range(len(self.seq)) if self.mat[i - 1][x]]
        if not scored_prev:
            return len(self.seq) + 1, len(self.seq) + 2
        return min(scored_prev) + self.min_len, max(scored_prev) + self.max_len

    def traceback(self, n_frag):
        score_li = []
        if not n_frag:
            for i in range(len(self.mat)):
                for j in range(len(self.seq)):
                    if self.mat[i][j] and j + self.max_len > len(self.seq):
                        score_li.append((self.mat[i][j], self.split_mat[i][j], i, j))
        else:
            i = n_frag - 1
            for j in range(len(self.seq)):
                if self.mat[i][j] and j + self.max_len > len(self.seq):
                    score_li.append((self.mat[i][j], self.split_mat[i][j], i, j))

        self.sorted_score_li = sorted(score_li, key=lambda x: -x[0])
        assert self.sorted_score_li, (
            'No starting point for traces. '
            'Usually due to infeasible fragment length and fragment number combination.'
        )

        self.traces, self.trace_scores, self.trace_split_scores = [], [], []

        for n, (score, split_scores, i, j) in enumerate(self.sorted_score_li[:self.n_trace], 1):
            trace = self.__trace(i, j)
            print('Trace', n, ' breaks:', ' ,'.join(str(x) for x in trace),
                  '; Total score:', score)
            self.traces.append(trace)
            self.trace_scores.append(score)
            self.trace_split_scores.append(split_scores)

    def __trace(self, i, j):
        js = [j]
        while i:
            key = self.make_key(i, j)
            i, j = self.unmake_key(self.trace_di[key])
            js.append(j)
        return js

    def write_outfile(self, filename_prefix):
        import json

        out_name = filename_prefix + '.txt'
        mat_name = filename_prefix + '_matrix.json'
        trace_name = filename_prefix + '_trace.json'

        outlines = []
        a = outlines.append

        a('File contains the formatted output of the generated fragments')

        with open(mat_name, 'w') as fp:
            json.dump(self.mat, fp)
        a('Matrix saved as: ' + mat_name)
        print('Wrote matrix file:', mat_name)

        with open(trace_name, 'w') as fp:
            json.dump(self.trace_di, fp)
        print('Wrote trace file:', trace_name)
        a('Trace saved as: ' + trace_name)
        a('')
        a('Input sequence')
        a('==============')
        a(self.seq)
        a('')
        a('Length of input: ' + str(len(self.seq)))

        a('Parameters Used')
        a('===============')
        a('\tLower limit of fragment length: {}'.format(self.min_len))
        a('\tUpper limit of fragment length: {}'.format(self.max_len))
        a('\tStart position index: {}'.format(self.start_site))
        a('\tEfficiency weight: {}'.format(self.eff_w))
        a('\tFidelity weight: {}'.format(self.fid_w))
        a('\tNumber of output traces: {}'.format(self.n_trace))
        a("\t5' overhang: {}".format(self._5overhang))
        a("\t3' overhang: {}".format(self._3overhang))
        a('\tInitiation score: {}'.format(self.inti_score))
        a('\tExcluded positions: {}'.format(' ,'.join(str(x) for x in self.exclusion_list)))
        a('\tAlien overhangs: {}'.format(' ,'.join(str(x) for x in self.alien_overhangs)))
        a('')
        a('Tracebacks (Predictions)')
        a(''.join(['-' for _ in range(len(self.seq))]))

        for i in range(len(self.traces)):
            a('Traceback ' + str(i + 1))
            a('Score (higher is better): {:3.2f}'.format(self.trace_scores[i]))
            a('Fragmentation positions (starts at 0): ' + ' ,'.join(str(x) for x in self.traces[i]))

            eff_tot, fid_tot = self.trace_split_scores[i]
            a('Net ligation efficiency: {:5.2f} %'.format(eff_tot * 100))
            a('Net Fidelity (accuracy): {:5.2f} %'.format(fid_tot * 100))
            a('Net success rate: {:5.2f} %'.format(eff_tot * fid_tot * 100))
            a('')

            a('Numbering 100s: ' + ''.join(
                [str(x)[-3] if (len(str(x)) >= 3 and str(x)[-1] == '0' and str(x)[-2] == '0') else ' '
                 for x in range(0, len(self.seq))]))
            a('Numbering 10s: ' + ''.join(
                [str(x)[-2] if (len(str(x)) >= 2 and str(x)[-1] == '0') else ' '
                 for x in range(0, len(self.seq))]))
            a('Numbering 1s: ' + ''.join([str(x)[-1] for x in range(0, len(self.seq))]))
            a("Sequence :5'-" + self.seq + "-3'")
            a('Break points: ' + ''.join(
                '^' if x in self.traces[i][:-1] else ' ' for x in range(len(self.seq))))

            fragments, effs, fids, range_li = self.get_fragments(self.traces[i])
            a('')
            a('Fragments')
            a('---------')
            for fi, fragment in enumerate(fragments):
                full_frag = self._5overhang + fragment + self._3overhang
                a('Fragment {}'.format(fi + 1))
                a("Range {}-{}; length: {}; length with additions: {}".format(
                    range_li[fi][0], range_li[fi][1], len(fragment), len(full_frag)))
                a("5' and 3' efficiencies: {}, {}; 5' and 3' fidelities: {}, {}".format(
                    effs[fi][0], effs[fi][1], fids[fi][0], fids[fi][1]))
                gap_spots = (
                    list(range(len(self._5overhang), len(self._5overhang) + 4))
                    + list(range(len(full_frag) - len(self._3overhang) - 4,
                                 len(full_frag) - len(self._3overhang)))
                )
                a('Overhang sites: ' + ''.join('m' if x in gap_spots else ' ' for x in range(len(full_frag))))
                a("Sequence :5'-{}-3'".format(full_frag))
                a("Reverse complement:3'-{}-5'".format(self.get_rc(full_frag, reverse=False)))
                a('Non ATCG bases: ' + ''.join('*' if x not in ['A', 'T', 'G', 'C'] else ' ' for x in full_frag))
                a('')

            a(''.join(['-' for _ in range(len(self.seq))]))

        with open(out_name, 'w') as fp:
            fp.writelines([x + '\n' for x in outlines])
        print('Wrote output file:', out_name)

    def get_fragments(self, js):
        js_rev = js[::-1]
        fragments = [self.seq[:js_rev[1] + 4]]
        eff, fid = self.get_score_for_j(js_rev[1])
        efficiencies = [('-', self.fl(eff))]
        fidelities = [('-', self.fl(fid))]
        range_li = [(0, js_rev[1])]

        for i in range(1, len(js_rev)):
            try:
                fragments.append(self.seq[js_rev[i]:js_rev[i + 1] + 4])
                eff, fid = self.get_score_for_j(js_rev[i])
                _eff, _fid = self.get_score_for_j(js_rev[i + 1])
                efficiencies.append((self.fl(eff), self.fl(_eff)))
                fidelities.append((self.fl(fid), self.fl(_fid)))
                range_li.append((js_rev[i], js_rev[i + 1]))
            except IndexError:
                break

        fragments.append(self.seq[js_rev[i]:])
        efficiencies.append((self.fl(_eff), '-'))
        fidelities.append((self.fl(_fid), '-'))
        range_li.append((js_rev[i], len(self.seq) + 1))

        return fragments, efficiencies, fidelities, range_li

    def get_score_for_j(self, j):
        ovrhg = self.seq[j:j + 4]
        return self.score_di[ovrhg]

    def get_rc(self, DNA, reverse=True):
        complements = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M',
            'R': 'Y', 'Y': 'R', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N'
        }
        if reverse:
            return ''.join(complements[x] for x in DNA[::-1])
        else:
            return ''.join(complements[x] for x in DNA)

    def make_key(self, i, j): return str(i) + '_' + str(j)
    def unmake_key(self, key): return int(key.split('_')[0]), int(key.split('_')[1])
    def round_up(self, number): return int(number) + (number % 1 > 0)
    def fl(self, a): return '{:5.2f}'.format(a)
