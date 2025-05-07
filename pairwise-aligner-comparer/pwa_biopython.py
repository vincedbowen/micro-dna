from Bio import Align

import argparse
from Bio import Align


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--A', required=True, help='First DNA sequence (query)')
    parser.add_argument('--B', required=True, help='Second DNA sequence (reference)')
    return parser.parse_args()


def align_sequences(seq_a: str, seq_b: str):
    aligner = Align.PairwiseAligner()
    # Match SSW-style scoring
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -1

    score = aligner.score(seq_a, seq_b)
    alignment = aligner.align(seq_a, seq_b)[0]

    return {
        "score": score,
        "alignment": alignment
    }


def main():
    args = get_args()
    result = align_sequences(args.A, args.B)

    print(f"Score: {result['score']}")
    print(result['alignment']) 


if __name__ == '__main__':
    main()
