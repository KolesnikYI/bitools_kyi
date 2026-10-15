import os

from helpers import constants, fastq


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: int | tuple[int, int] = constants.DEFAULT_GC_BOUNDS,
    # length_bounds: tuple[int, int] = constants.DEFAULT_LENGTH_BOUNDS,
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
            fastq.is_in_gc_bounds(
                sequence=one_read[1],
                gc_bounds=fastq.to_tuple(gc_bounds),
            )

            [output_file.write(line) for line in one_read]


if __name__ == "__main__":
    filter_fastq("test_input_fastq.fastq", "test_output.fastq", gc_bounds=40)
