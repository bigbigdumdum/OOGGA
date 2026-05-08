"""
Command-line interface for OOGGA.

Fragments a DNA sequence optimally for Golden Gate assembly.

Copyright (C) 2025 Mukundan S <mukundan.kollam@gmail.com>
Licensed under CC BY-NC-SA 4.0.
Please cite https://doi.org/10.1101/2025.06.16.659877
"""

import argparse
import sys
from .core import read_fasta, Dyna_frag


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Program takes in a DNA sequence and fragments it for Golden Gate cloning.
Overhangs of each fragment are chosen to minimise error rate and maximise ligation efficiency.
Scoring data from Potapov et al. 2018 (https://doi.org/10.1021/acssynbio.8b00333).

Helpful BsaI attachment sequences:
  3' attachment: GGCTACGGTCTCC
  5' attachment: GGAGACCGTAGCC
""",
    )

    parser.add_argument("input", type=str,
                        help="Input FASTA file (only the first sequence is processed)")
    parser.add_argument("min_len", type=int,
                        help="Minimum fragment length")
    parser.add_argument("max_len", type=int,
                        help="Maximum fragment length")
    parser.add_argument("output", type=str,
                        help="Output file prefix (directory will be created if needed)")

    parser.add_argument("-n_frag", type=int, default=False,
                        help="Fixed number of fragments (omit for auto-detection)")
    parser.add_argument("-data_file", type=str,
                        default="./lib/FileS04_T4_18h_37C.csv",
                        help="Path to Potapov et al. 2018 CSV scoring table. "
                             "Download FileS04_T4_18h_37C.xlsx from "
                             "https://pubs.acs.org/doi/suppl/10.1021/acssynbio.8b00333/suppl_file/sb8b00333_si_002.zip "
                             "and convert to CSV.")
    parser.add_argument("-start", type=int, default=0,
                        help="Start position index of the first fragment")
    parser.add_argument("-eff_w", type=float, default=1,
                        help="Exponent weight for ligation efficiency metric")
    parser.add_argument("-fid_w", type=float, default=1,
                        help="Exponent weight for fidelity metric")
    parser.add_argument("-ovrhg_iden", type=int, default=2,
                        help="Maximum allowed identity between two overhangs")
    parser.add_argument("-n_sol", type=int, default=5,
                        help="Number of solutions to output")
    parser.add_argument("-exclude", nargs='*', default=[], type=int,
                        help="DNA indices where overhangs cannot be placed")
    parser.add_argument("-add3", type=str, default='',
                        help="Sequence to add to the 3' end (e.g. restriction enzyme site)")
    parser.add_argument("-add5", type=str, default='',
                        help="Sequence to add to the 5' end (e.g. restriction enzyme site)")
    parser.add_argument("-alien_overhangs", nargs='*', default=[], type=str,
                        help="External overhangs to include in identity checking")

    a = parser.parse_args()

    if a.data_file is None:
        print(
            "Confused OOGGA booga: -data_file is required.\n\n"
            "Download the scoring table from Potapov et al. 2018:\n"
            "  https://pubs.acs.org/doi/suppl/10.1021/acssynbio.8b00333/suppl_file/sb8b00333_si_002.zip\n\n"
            "Extract the zip, open FileS04_T4_18h_37C.xlsx in Excel or LibreOffice,\n"
            "and save it as FileS04_T4_18h_37C.csv. Then pass it with:\n"
            "  oogga ... -data_file /path/to/FileS04_T4_18h_37C.csv\n",
            file=sys.stderr,
        )
        sys.exit(1)

    assert (a.max_len - a.min_len) >= 4, 'Confused OOGGA booga: Difference between minimum and maximum length must be greater than 4.'

    import os
    out_dir = os.path.dirname(a.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    try:
        seq = list(read_fasta(a.input).values())[0]
        df = Dyna_frag(
            seq, a.min_len, a.max_len, a.data_file,
            n_frag=a.n_frag, n_trace=a.n_sol,
            start_site=a.start, eff_w=a.eff_w, fid_w=a.fid_w,
            exclusion_list=a.exclude, inti_score=1,
            max_overhang_identity=a.ovrhg_iden,
            alien_overhangs=a.alien_overhangs,
        )
        df._5overhang = a.add3
        df._3overhang = a.add5
        df.write_outfile(a.output)
    except Exception as e:
        print('Confused OOGGA booga:', e)
        sys.exit(1)


if __name__ == '__main__':
    main()
