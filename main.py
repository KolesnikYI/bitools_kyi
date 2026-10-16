import os

from helpers import constants, fastq


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: int | tuple[int, int] = constants.DEFAULT_GC_BOUNDS,
    length_bounds: tuple[int, int] = constants.DEFAULT_LENGTH_BOUNDS,
    quality_threshold: int = constants.DEFAULT_QUALITY_THRESHOLD,
) -> None:
    """Function for filtering fastq file by quality settings.

    Args:
        input_fastq (str): input file path.
        output_fastq (str): output file name.
        gc_bounds (tuple[int, int], optional): tuple of gc bounds, also may set the one number - it will be upper bound. Defaults to (0, 100).
        length_bounds (tuple[int, int], optional): tuple of bounds length, also may set the one number - it will be max length size. Defaults to (0, 2**32).
        quality_threshold (int, optional): quality threshold (phred33 scale). Defaults to 0.
    """
    output_path = constants.FILTERED_DIR + os.sep + output_fastq
    with open(input_fastq, "r") as f, open(output_path, "w") as output_file:
        if not os.path.exists(constants.FILTERED_DIR):
            os.makedirs(constants.FILTERED_DIR)

        while True:
            header, sequence, plus, quality = [f.readline() for _ in range(4)]
            if not header:
                break
            one_read = (header, sequence, plus, quality)
            if fastq.keep_in_file(
                one_read,
                gc_bounds,
                length_bounds,
                quality_threshold,
            ):
                [output_file.write(line) for line in one_read]


def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: str = "",
) -> None:
    """Convert multiple line fasta file to monosequence line file.

    Args:
        input_fasta (str): path to fasta file
        output_fasta (str, optional): name of output file, default will "convertered_input_file".
    """
    if not output_fasta:
        output_fasta = "convertered_" + input_fasta
    with open(input_fasta, "r") as f, open(output_fasta, "w") as output_file:
        chunk = []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if chunk:
                    output_file.write("".join(chunk) + "\n")
                    chunk.clear()
                output_file.write(line + "\n")
            else:
                chunk.append(line)

        if chunk:
            output_file.write("".join(chunk) + "\n")


def parse_blast_output(
    input_file: str,
    output_file: str,
) -> None:
    """Parsing proteins description from blast output.

    Args:
        input_file (str): path to blast_output file for parsing
        output_file (str): name of ouptut parsed file
    """
    with open(input_file, "r") as f, open(output_file, "w") as output:
        all_lines = f.readlines()
        result = []
        for index, line in enumerate(all_lines, start=1):
            if constants.START_STRING_PARSING in line:
                target_line = all_lines[index + 2]
                result.append(
                    target_line[: constants.END_LINE_INDEX].strip() + "\n"
                )
        result.sort()
        output.write("".join(result) + "\n")
