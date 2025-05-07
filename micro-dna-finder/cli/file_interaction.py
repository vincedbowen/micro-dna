import pysam
import csv
from beautiful_printer import status_message

def load_bam(bam_path: str) -> pysam.AlignmentFile:
    """Loads BAM file specified by user into memory

    Args:
        bam_path (str): path and filename to be loaded

    Returns:
        pysam.AlignmentFile: processed bam Alignment File
    """
    return pysam.AlignmentFile(bam_path, "rb")

def load_fasta(fasta_path: str) -> str:
    """Loads Fasta file specified by user into memory

    Args:
        fasta_path (str): path and filename to be loaded

    Returns:
        fasta.fetch(contig): processed fasta Reference file. In this case it is just the first chromosome
    """
    fasta = pysam.FastaFile(fasta_path)
    contig = fasta.references[0]
    return fasta.fetch(contig)

@status_message("Saving results to txt file")
def save_results_txt(results: list, filename: str) -> None:
    """Writes results to a txt file in a human-readable format

    Args:
        results (list): list of dictionaries containing high evidence of a circular DNA sequence
        filename (str): path and filename to be saved (txt format)
    """
    with open(filename, "w") as f:
        for r in results:
            f.write("COMPARING SOFT CLIPS\n")
            f.write(f"Start: {r['start_pos']} End: {r['end_pos']}\n")
            f.write(f"Score: {r['read_read_score']}\n")
            f.write(str(r['read_read_alignment']) + "\n")

            f.write("COMPARING START WITH REFERENCE\n")
            f.write(f"Score: {r['start_ref_score']}\n")
            f.write(str(r['start_ref_alignment']) + "\n")

            f.write("COMPARING END WITH REFERENCE\n")
            f.write(f"Score: {r['end_ref_score']}\n")
            f.write(str(r['end_ref_alignment']) + "\n")
            f.write(f"Evidence Score: {r['evidence_score']}\n")
            f.write("-" * 75 + "\n\n")
            
@status_message("Saving results to tsv/csv file")
def save_results_sv(results: list, filename: str) -> None:
    """Writes results to a tsv file in a human-readable format, but more optimized for machine readability.
    Doesn't include the alignment information, but includes the scores and positions.

    Args:
        results (list): list of dictionaries containing high evidence of a circular DNA sequence
        filename (str): path and filename to be saved (tsv or csv format)
    """
    fieldnames = [
        "start_pos",
        "end_pos",
        "read_read_score",
        "start_ref_score",
        "end_ref_score",
        "evidence_score",
    ]
    # Write the results to a CSV or TSV file
    if filename.endswith('.csv'):
        with open(filename, "w", newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow({field: r.get(field, "") for field in fieldnames})
    else:
        with open(filename, "w", newline='') as tsvfile:
            writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for r in results:
                writer.writerow({field: r.get(field, "") for field in fieldnames})

def clear_file(filename: str) -> None:
    """Clears the contents of a file

    Args:
        filename (str): path and filename to be cleared
    """
    open(filename, "w").close()