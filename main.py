import os

from helpers import constants


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: int | tuple[int, int] = constants.DEFAULT_GC_BOUNDS,
    length_bounds: tuple[int, int] = constants.DEFAULT_LENGTH_BOUNDS,
    # quality_threshold: int = constants.DEFAULT_QUALITY_THRESHOLD,
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
            is_in_gc_bounds(
                sequence=one_read[1],
                gc_bounds=to_tuple(gc_bounds),
            )

            [output_file.write(line) for line in one_read]


def is_in_length_bounds(
    sequence: str,
    length_bounds: tuple[int, int],
) -> bool:
    """Check that sequence length in filter parameter.

    Args:
        sequence (str): line with sequence from one read.
        length_bounds (tuple[int, int]): filter bounds.

    Returns:
        bool: true or false
    """

    return (
        length_bounds[constants.LOWER_BOUND_INDEX]
        <= len(sequence)
        < length_bounds[constants.UPPER_BOUND_INDEX]
    )


def to_tuple(input_value: int | tuple[int, int]) -> tuple[int, int]:
    """Format int to (0, int) format.

    Args:
        input_value (int | tuple[int, int]): may be tuple, or one number.

    Returns:
        tuple[int, int]: correct tuple
    """
    if isinstance(input_value, int):
        return (0, input_value)
    return input_value


def is_in_gc_bounds(
    sequence: str,
    gc_bounds: tuple[int, int],
) -> bool:
    """Check GC fraction in filter parameter.

    Args:
        sequence (str): line with sequence from one read.
        gc_bounds (tuple[int, int]): filter bounds.

    Returns:
        bool: true or false
    """
    gc_fraction = (
        (sequence.count("G") + sequence.count("C")) / len(sequence)
    ) * 100
    return (
        gc_bounds[constants.LOWER_BOUND_INDEX]
        < gc_fraction
        < gc_bounds[constants.UPPER_BOUND_INDEX]
    )


if __name__ == "__main__":
    filter_fastq("test_input_fastq.fastq", "test_output.fastq", gc_bounds=40)
