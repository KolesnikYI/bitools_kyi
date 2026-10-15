from helpers import constants


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
