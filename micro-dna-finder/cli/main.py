from utils import *
from file_interaction import *

import argparse
from colorama import Fore

def get_args():
    """Handles CLI arguments

    Returns:
        args: The CLI arguments as key-value pairs
    """    
    parser = argparse.ArgumentParser()
    parser.add_argument('-bf',
                        '--bam_file',
                        default=None,
                        type=str,
                        required=True,
                        help='BAM file')
    parser.add_argument('-ff',
                        '--fasta_file',
                        default=None,
                        type=str,
                        required=True,
                        help='Fasta file')
    parser.add_argument('-of',
                        '--output_file',
                        default=None,
                        type=str,
                        required=True,
                        help='Output file in txt, tsv, or csv format')
    return parser.parse_args()

def main():
    args = get_args()
    # Load in files for analysis
    bam = load_bam(args.bam_file)
    contig = load_fasta(args.fasta_file)
    clear_file(args.output_file)
    # Locate start soft-clips and filter high support sites
    start_soft_clips = find_soft_clips(bam, start=True)
    start_hsc = find_high_support_sites(start_soft_clips)
    # Locate end soft-clips and filter high support sites
    end_soft_clips = find_soft_clips(bam, start=False)
    end_hsc = find_high_support_sites(end_soft_clips)
    # Align soft-clips with each other
    aligned_sc = align_soft_clips(start_hsc, end_hsc)
    # Generate final results by aligning soft-clips with contig and sorting
    results = generate_final_results(aligned_sc, contig)
    if args.output_file.endswith('.tsv') or args.output_file.endswith('.csv'):
        save_results_sv(results, args.output_file)
    elif args.output_file.endswith('.txt'):
        save_results_txt(results, args.output_file)

if __name__ == "__main__":
    main()

