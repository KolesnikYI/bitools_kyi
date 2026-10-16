from helpers import constants
from statistics import mean


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


def is_above_quality_threshold(
    quality: str,
    threshold: int,
) -> bool:
    """Check quolity line for threshold filter

    Args:
        quality (str): line with quolity encryption from one read.
        threshold (int): threshold value

    Returns:
        bool: true or false
    """
    int_read_quality = []
    for symbol in quality:
        int_read_quality.append(ord(symbol) - constants.PHRED_CONSTANT)
    return mean(int_read_quality) >= threshold


def keep_in_file(
    one_read: tuple[str, ...],
    gc_bounds: int | tuple[int, int],
    length_bounds: int | tuple[int, int],
    quality_threshold: int,
) -> bool:
    """Union all filters in one place

    Args:
        one_read (tuple[str, ...]): tuple with read components.
        gc_bounds (int | tuple[int, int]): bounds of GC fraction.
        length_bounds (int | tuple[int, int]): length of sequence.
        quality_threshold (int): quality threshold.

    Returns:
        bool: true if all filter passed
    """
    sequence = one_read[constants.SEQUENCE_INDEX]
    quality = one_read[constants.QUALITY_INDEX]
    gc_bounds = to_tuple(gc_bounds)
    length_bounds = to_tuple(length_bounds)
    return (
        is_in_gc_bounds(sequence, gc_bounds)
        and is_in_length_bounds(sequence, length_bounds)
        and is_above_quality_threshold(quality, quality_threshold)
    )
