"""
Command-line interface for eval_frags - evaluate overhang scoring.

Evaluates a list of 4-bp overhangs using the OOGGA scoring scheme.

Copyright (C) 2025 Mukundan S <mukundan.kollam@gmail.com>
Licensed under CC BY-NC-SA 4.0.
Please cite https://doi.org/10.1101/2025.06.16.659877
"""

import argparse
import sys
from .core import read_fasta, load_csv_table_as_di


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Evaluates a list of 4-bp overhang sequences and reports
net ligation efficiency and probability of correct assembly.
Scoring data from Potapov et al. 2018 (https://doi.org/10.1021/acssynbio.8b00333).
""",
    )

    parser.add_argument("overhangs", nargs='*', type=str,
                        help="4-bp overhang sequences (space-separated, 5'->3')")
    parser.add_argument("-data_file", type=str, default=None,
                        help="Path to Potapov et al. 2018 CSV scoring table.")
    parser.add_argument("-fasta_input", type=str, default=False,
                        help="FASTA input (e.g. SplitSet output) to read overhangs from")

    a = parser.parse_args()

    if a.data_file is None:
        print(
            "Confused OOGGA booga: -data_file is required.\n\n"
            "Download from:\n"
            "  https://pubs.acs.org/doi/suppl/10.1021/acssynbio.8b00333/suppl_file/sb8b00333_si_002.zip\n"
            "Convert FileS04_T4_18h_37C.xlsx to CSV and pass with -data_file.\n",
            file=sys.stderr,
        )
        sys.exit(1)

    ovrhgs = list(a.overhangs)

    if a.fasta_input:
        seq_di = read_fasta(a.fasta_input)
        for x in seq_di:
            ovrhgs.append(x[-4:])

    print('Fragments:', ', '.join(str(x) for x in ovrhgs))

    score_di = load_csv_table_as_di(a.data_file)
    net_eff, net_fid = 1, 1

    print('{:10s} {:12s} {:10s}'.format('Overhang', 'Efficiency(%)', 'Fidelity(%)'))
    for ovrhg in ovrhgs:
        try:
            eff, fid = score_di[ovrhg]
        except KeyError:
            print('{:10s}  NOT FOUND in scoring table'.format(ovrhg))
            continue
        print('{:10s} {:10.2f} {:10.2f}'.format(ovrhg, round(eff), round(fid)))
        net_eff = net_eff * (eff / 100)
        net_fid = net_fid * (fid / 100)

    print('\nNet efficiency = {:.4f} %'.format(net_eff * 100))
    print('Net fidelity   = {:.2f} %'.format(net_fid * 100))


if __name__ == '__main__':
    main()
