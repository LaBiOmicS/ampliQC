"""
FASTQ streaming parser wrapper using high performance dnaio / xopen.
"""

from pathlib import Path
from typing import Generator, Tuple, Union, Optional
import dnaio


class FastqRecord:
    __slots__ = ("name", "sequence", "qualities")

    def __init__(self, name: str, sequence: str, qualities: str):
        self.name = name
        self.sequence = sequence
        self.qualities = qualities


def read_fastq(
    filepath: Union[str, Path]
) -> Generator[FastqRecord, None, None]:
    """
    Stream FASTQ records from uncompressed or gzipped files using dnaio.
    """
    filepath = str(filepath)
    with dnaio.open(filepath, mode="r") as f:
        for record in f:
            yield FastqRecord(
                name=record.name,
                sequence=record.sequence,
                qualities=record.qualities,
            )


def read_fastq_paired(
    r1_path: Union[str, Path],
    r2_path: Union[str, Path]
) -> Generator[Tuple[FastqRecord, FastqRecord], None, None]:
    """
    Stream paired-end FASTQ records from R1 and R2.
    """
    with dnaio.open(str(r1_path), file2=str(r2_path), mode="r") as f:
        for r1, r2 in f:
            yield (
                FastqRecord(r1.name, r1.sequence, r1.qualities),
                FastqRecord(r2.name, r2.sequence, r2.qualities),
            )
