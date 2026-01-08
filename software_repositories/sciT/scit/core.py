import pysam
import numpy as np
from typing import Literal

FASTQ_MODALITY_DATATYPES = Literal['RNA']
MULTIPLEXER_ALIASES = Literal['sciT']


def wasp_read_is_assigned_to_haplotype_1(read: pysam.AlignedSegment) -> bool:
    return read.has_tag('vW') and \
        read.get_tag('vW') == 1 and \
        read.has_tag('vA') and \
        np.median( read.get_tag('vA') ).astype(int) == 1

def wasp_read_is_assigned_to_haplotype_2(read: pysam.AlignedSegment) -> bool:
    return read.has_tag('vW') and \
        read.get_tag('vW') == 1 and \
        read.has_tag('vA') and \
        np.median( read.get_tag('vA') ).astype(int) == 2
