import argparse
from skbio.alignment import StripedSmithWaterman


def get_args():
    parser = argparse.ArgumentParser(description="Smith-Waterman alignment using skbio's StripedSmithWaterman.")
    parser.add_argument('--A', required=True, help='First DNA sequence (query)')
    parser.add_argument('--B', required=True, help='Second DNA sequence (target)')
    return parser.parse_args()


def align_sequences(seq_a: str, seq_b: str):
    aligner = StripedSmithWaterman(seq_a) 
    result = aligner(seq_b)

    return {
        "score": result.optimal_alignment_score,
        "start_query": result.query_begin,
        "end_query": result.query_end,
        "cigar": result.cigar,
        "alignment": result.aligned_query_sequence,
        "reference": result.aligned_target_sequence
    }


def main():
    args = get_args()
    result = align_sequences(args.A, args.B)

    print(f"Score: {result['score']}")
    print(f"Query  : {result['alignment']}")
    print(f"Target : {result['reference']}")
    print(f"CIGAR  : {result['cigar']}")


if __name__ == '__main__':
    main()
