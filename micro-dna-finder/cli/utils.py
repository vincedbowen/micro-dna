import pysam
from Bio import Align
from alignment_models import AlignmentDetails, AlignerDetails
from colorama import Fore
import logging
from beautiful_printer import status_message
logger = logging.getLogger(__name__)

@status_message("Locating Soft-Clips")
def find_soft_clips(bam: pysam.AlignmentFile, start: bool) -> dict[int, AlignmentDetails]:
    """Finds soft-clipped reads in a BAM file

    Args:
        bam (pysam.AlignmentFile): the alignment file to find soft-clips in
        start (bool): flag to indicate if the soft-clip is at the start or end of the read

    Returns:
        dict[int, AlignmentDetails]: 1-based position of the soft-clip in the reference genome and the corresponding sequence, support count, and length of the soft-clip for start soft-clips
    """
    # In the format {position: [clipped_seq]}
    scr: dict[int, AlignmentDetails] = {} 
    # Iterate through each read in the BAM file, resetting the iterator every time
    for read in bam.fetch():
        cigars = read.cigartuples
        seq = read.query_sequence
        pos = read.reference_start + 1
        # Soft-clip at start
        if cigars[0][0] == 4 and cigars[-1][0] == 0 and start:
            if pos not in scr:
                scr[pos] = AlignmentDetails(seq, sc_len=cigars[0][1])
            else:
                scr[pos].update_count()   
        # Soft-clip at end
        elif cigars[0][0] == 0 and cigars[-1][0] == 4 and not start:
            if pos not in scr:
                scr[pos] = AlignmentDetails(seq)
            else:
                scr[pos].update_count()
    return scr

@status_message("Filtering high-support soft-clips")
def find_high_support_sites(soft_clips_dict: dict[int, AlignmentDetails], min_support: int = 10) -> dict[int, AlignmentDetails]:
    """Finds high-support soft-clip sites in a dictionary of soft-clips

    Args:
        soft_clips_dict (dict[int, AlignmentDetails]): the dictionary of soft-clips
        min_support (int, optional): Arbitrary number defined as the threshold for high support. Defaults to 10.

    Returns:
        dict[int, AlignmentDetails]: the filtered dictionary of soft-clips with high support
    """
    hi_supp: dict[int, AlignmentDetails] = {} 
    for key, val in soft_clips_dict.items():
        if val.count >= min_support:
            hi_supp[key] = val
    return hi_supp

def get_seq_from_pos(start_pos: int, end_post: int, contig: str) -> str:
    """Gets the sequence from a contig at specified positions

    Args:
        start_pos (int): start position of the sequence
        end_post (int): end position of the sequence
        contig (str): the processed contig

    Returns:
        str: the sequence from the contig at the specified positions
    """
    seq = contig[start_pos:end_post]
    return seq

def aligner_init(match_score: int = 2, mismatch_score: int = -1, open_gap_score: int = -0.5, extend_gap_score: int = -0.1, mode: str = 'global', sc: bool = False) -> Align.PairwiseAligner:
    """Initializes the pairwise aligner with the specified parameters

    Args:
        match_score (int, optional): match score defined in BioPython. Defaults to 2.
        mismatch_score (int, optional): mismatch score defined in BioPython. Defaults to -1.
        open_gap_score (int, optional): open gap score defined in BioPython. Defaults to -0.5.
        extend_gap_score (int, optional): extend gap score defined in BioPython. Defaults to -0.1.
        mode (str, optional): 'global' or 'local' mode which determines the aligner methods. Defaults to 'global'.
        sc (bool, optional): flag to determine if the alignment is on two soft-clips. Defaults to False.

    Returns:
        Align.PairwiseAligner: the initialized pairwise aligner
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    if sc:
        aligner.query_left_open_gap_score = -0.1
        aligner.query_right_open_gap_score = -0.1
        aligner.query_left_extend_gap_score = -0.1
        aligner.query_right_extend_gap_score = -0.1 
    # logger.info(f"Pairwise sequence alignment initialized using the {Fore.BLUE}{aligner.algorithm}{Fore.RESET} with {Fore.BLUE}{aligner.mode}{Fore.RESET} mode.")
    return aligner

def align_strs(aligner: Align.PairwiseAligner, query: str, target: str) -> AlignerDetails:
    """Aligns two sequences in string format using the specified aligner

    Args:
        aligner (Align.PairwiseAligner): the pairwise aligner
        query (str): the query sequence
        target (str): the target sequence (another soft-clip or the reference genome)

    Returns:
        AlignerDetails: the alignment details including the alignment sequence and score
    """
    alignments = aligner.align(query, target)
    # If no alignments are found, return None and a score of -inf to make evidence score calculation alter
    if not alignments:
        return AlignerDetails(None, float("-inf"))
    alignment = alignments[0]
    score = alignment.score
    return AlignerDetails(alignment, score)

@status_message("Aligning starting soft-clips with ending soft-clips")
def align_soft_clips(start_hsc: dict, end_hsc: dict) -> list[dict]:
    """Brute-force aligns all starting soft-clips with all ending soft-clips

    Args:
        start_hsc (dict): dictionary of starting soft-clips
        end_hsc (dict): dictionary of ending soft-clips

    Returns:
        list[dict]: list of dictionaries containing the alignment details for each pair of starting and ending soft-clips
    """
    results = []
    a = aligner_init(2, -2, -1, -0.5, sc=True)
    for start_pos, start_value in start_hsc.items():
        start_seq = start_value.seq
        for end_pos, end_value in end_hsc.items():
            end_seq = end_value.seq
            # Get best alignment
            res = align_strs(a, start_seq, end_seq)
            results.append({
                "start_pos": start_pos,
                "end_pos": end_pos,
                "score": res.score,
                "alignment": res.alignment,
                "start_seq":start_seq,
                "end_seq":end_seq,
                "start_len": start_value.sc_len,
            })
    
    return results

def align_with_reference(reference: str, position: int, short_clip: str, aligner: Align.PairwiseAligner, start_len: int = 0) -> dict:
    """Aligns a soft-clip with the reference genome at a specified position

    Args:
        reference (str): the reference genome sequence
        position (int): the position of the soft-clip in the reference genome
        short_clip (str): the soft-clip sequence
        aligner (Align.PairwiseAligner): the pairwise aligner
        start_len (int, optional): the length of the soft-clip if it is a start soft-clip. Defaults to 0.

    Returns:
        dict: a dictionary containing the alignment details including the alignment sequence and score
    """
    position = position - 1 - start_len
    ref_start_seq = reference[position:min(position + 42, len(reference) - 1)]
    # Align soft-clip with reference start sequence
    start_alignment = align_strs(aligner, ref_start_seq, short_clip)
    # Align soft-clip with reference end sequence
    result = {
        "alignment": str(start_alignment.alignment),
        "score": start_alignment.score,
        "reference_sequence": ref_start_seq
    }
    return result

def create_evidence_score(sc_alignment: dict, start_ref_alignment: dict, end_ref_alignment: dict, ideal: int = 200) -> float:
    """Creates a score for the evidence based on the alignment scores of the soft-clip and the reference genome using weighted scoring and a quadratic penalty for the span of the soft-clip.
    The score is calculated as follows:
    score = (short_clip_score * 8 + start_ref_score * 2 + end_ref_score * 2) / penalty
    where penalty is a quadratic function of the span of the soft-clip.
    The penalty is calculated as:
    penalty = ((span - ideal) / ideal) ** 2 + 1
    where span is the length of the soft-clip and ideal is the (arbitrary) ideal length of the soft-clip.

    Args:
        sc_alignment (dict): dictionary containing the alignment details of the soft-clip
        start_ref_alignment (dict): the alignment details of the start soft-clip and reference genome
        end_ref_alignment (dict): the alignment details of the end soft-clip and reference genome
        ideal (int, optional): the ideal length of the soft-clip. Defaults to 200.

    Returns:
        float: _description_
    """
    span = sc_alignment["end_pos"] - sc_alignment["start_pos"]
    # Quadratic function where the minimum penalty is at the ideal length and the penalty increases as the span deviates from the ideal length
    penalty = ((span - ideal) / ideal) ** 2 + 1
    score = (
        sc_alignment["score"] * 8 +
        start_ref_alignment["score"] * 2 +
        end_ref_alignment["score"] * 2
    ) / penalty
    return score

@status_message("Aligning soft-clips with reference genome")
def generate_final_results(aligned_sc:list[dict], contig: str) -> list[dict]:
    """Aligns soft-clips with the reference genome and filters out low-scoring alignments
    and short soft-clips. The alignments are sorted by the evidence score in descending order.

    Args:
        aligned_sc (list[dict]): A list of dictionaries containing the alignment details of the soft-clips
        contig (str): The reference genome sequence

    Returns:
        list[dict]: the filtered and sorted list of dictionaries containing the final alignment details
    """
    results = []
    # Initialize the aligner for aligning soft-clips with the reference genome
    aligner = aligner_init(2, -3, -25, -6, 'global')
    for aligned in aligned_sc:
        # Filter out low-scoring alignments and short soft-clips
        if aligned['score'] > 50 and aligned['end_pos'] - aligned['start_pos'] > 5:
            # Align the start and end soft-clips with the reference genome
            ssc = align_with_reference(
                reference=contig,
                position=aligned['start_pos'],
                short_clip=aligned['start_seq'],
                aligner=aligner,
                start_len = aligned['start_len']
            )
            esc = align_with_reference(
                reference=contig,
                position=aligned['end_pos'],
                short_clip=aligned['end_seq'],
                aligner=aligner
            )
            results.append({
                "start_pos": aligned['start_pos'],
                "end_pos": aligned['end_pos'],
                "read_read_score": aligned['score'],
                "read_read_alignment": str(aligned['alignment']),
                "start_ref_score": ssc['score'],
                "start_ref_alignment": ssc['alignment'],
                "end_ref_score": esc['score'],
                "end_ref_alignment": esc['alignment'],
                "evidence_score": create_evidence_score(aligned, ssc, esc)
            })
    # Sort the results by the evidence score in descending order
    results = sorted(results, key=lambda x: x["evidence_score"], reverse=True)
    return results

