


import subprocess
from typing import Callable, Dict, List, Tuple, Set
from collections import defaultdict, Counter
import pysam
from singlecellmultiomics.bamProcessing.bamFunctions import get_read_group_to_sample_dict
import argparse
from .transcriptome_bam_to_loom import generic_single_end_tx_read_pass_function


def count_unique_molecules_per_sample(bam_path: str, pass_function: Callable[[pysam.AlignedSegment], bool]) -> Tuple[Counter[str], Dict[str, str]] :
    """
    Returns the number of reads for each read group, and a read group dict: (read_group -> sample)
    { 'RNA.AAGJMJ3M5.1.Even2Bo2_Odd2Bo49.MB-10K_multiome_144': 'MB-10K_multiome_144', ... }
    """
    molecules: Counter[str] = Counter()
    with pysam.AlignmentFile(bam_path, 'rb', threads=3) as aln:
        read_groups = get_read_group_to_sample_dict(aln) #( RG -> sample)
        for read in aln:
            if not pass_function(read):
                continue
            sample = read_groups.get(read.get_tag('RG'))
            if sample is None: 
                continue
            molecules[sample] += 1
    return molecules, read_groups


def sciT_cell_filter(
            transcriptome_input_bam: str,
            transcriptome_output_bam: str,
            cell_min_transcriptome_count: int = 0) -> Set[str]:
    """
    Filter cells based on the number of molecules in the cell.
    When a cell passes any of the thresholds the cell is kept. 
    Returns a set with all samples which are kept
    """    
    cells_to_keep: set[str] = set()
    transcriptome_molecules, transcriptome_read_groups = count_unique_molecules_per_sample(transcriptome_input_bam, generic_single_end_tx_read_pass_function)
    for cell, count in transcriptome_molecules.items():
        print(cell, count, 'tx')
        if count >= cell_min_transcriptome_count:
            cells_to_keep.add(cell)
    subset_samples_from_bam(transcriptome_input_bam, transcriptome_output_bam, cells_to_keep, transcriptome_read_groups)
    return cells_to_keep

# Invert the readgroups to get SAMPLE-> RGs
def sample_to_read_group_ids(read_groups: Dict[str,str]) -> Dict[str,List[str]]:
    sample_to_rg: Dict[str,List[str]] = defaultdict(list)
    for rg, sample in read_groups.items():
        sample_to_rg[sample].append(rg)
    return sample_to_rg

def write_rg_file(path: str, keep_samples: Set[str], sample_to_rg: Dict[str,List[str]]):
    with open(path,'w') as o:
        for sample in keep_samples:
            for rg in sample_to_rg[sample]:
                o.write(rg+'\n')
def write_sample_file(path: str, keep_samples: Set[str]):
    with open(path,'w') as o:
        for sample in keep_samples:
            o.write(sample+'\n')
                                
                
def subset_samples_from_bam(in_bam: str, out_bam:str, keep_samples: Set[str], read_groups: Dict[str,str], n_compression_threads:int=4):
    sample_to_rg = sample_to_read_group_ids(read_groups)
    rgpath = out_bam.replace('.bam','.rg.txt')
    samplelist_path = out_bam.replace('.bam','.samples.txt')
    assert rgpath!=out_bam
    assert samplelist_path != out_bam
    write_sample_file(samplelist_path, keep_samples)
    write_rg_file(rgpath, keep_samples, sample_to_rg)
    cmd = ['samtools', 'view', '-@', str(n_compression_threads), '-o', out_bam, '-R', rgpath, in_bam, '--write-index']
    subprocess.run(cmd)
    
    
def run():
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="sciT low count filter"
    )
    argparser.add_argument(
        '--transcriptome_input_bam',
        type=str, help='Transcriptome bam file'
        )
    argparser.add_argument(
        '--transcriptome_output_bam',
        type=str, help='Transcriptome bam file'
        )
    argparser.add_argument(
        '--cell_min_transcriptome_count',
        default=10,
        type=int, help='Minimum number of transcriptome molecules per cell (qc passing and unique)'
        )
    
    args = argparser.parse_args()

    print(len(sciT_cell_filter(
            transcriptome_input_bam = args.transcriptome_input_bam,
            transcriptome_output_bam = args.transcriptome_output_bam,
            cell_min_transcriptome_count = args.cell_min_transcriptome_count)), 'cells passed the thresholds')
        
        
    