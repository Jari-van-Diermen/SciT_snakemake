#!/usr/bin/env python3
import argparse
import xopen
from typing import Iterator, Tuple,List
import yaml
from .core import FASTQ_MODALITY_DATATYPES, MULTIPLEXER_ALIASES
from .sciT_umi import autodetect_umi_len, UmiExtractor
import pysam

__version__ = '8 '

def extract_barcodes(reads_in: List[Tuple[pysam.FastxRecord, pysam.FastxRecord]], prefix:str, umi_extractor:UmiExtractor,  write_ci: bool=False):
    header = reads_in[0].name
    instrument, run_number, flow_cell_id, lane, tile, cluster_x_pos, cluster_y_pos, _, barcodeIdentification_header = header.replace(' ',':').split(':')
    dpm, y, even, odd = tuple(barcodeIdentification_header[1:-1].split(']['))
    # Calculate a cell-index
    even = even if even!='NOT_FOUND' else None
    odd = odd if odd != 'NOT_FOUND' else None
    y = y if y != 'NOT_FOUND' else None
    
    if even is not None and odd is not None and y is not None:
        # Calculate barcode index value:
        ntotal_barcodes = 96 # When this value is exceeded the script raises an error, an overestimate is always safe
        even_index = int(even.split('Bo')[-1])-1
        odd_index = int(odd.split('Bo')[-1])-1
        y_index =  int(y.split('Bot')[-1].split('_')[0])-1 # format is NYBot[0-9]_Stg
        assert even_index < ntotal_barcodes and odd_index < ntotal_barcodes and y_index < ntotal_barcodes, f'Barcode index out of range: {even_index}, {odd_index}, {y_index}'
        barcode_index = even_index*ntotal_barcodes + odd_index + (ntotal_barcodes**2)*y_index
        barcode_string = f'BC:{even}_{odd}_{y};bi:{barcode_index};'
        if write_ci:
            barcode_string += f'ci:{barcode_index};'
    else:
        # Compose the barcode as the parts which are set:
        barcode_string = f'bc:{"_".join([bc for bc in [even,odd,y] if bc is not None])};'

    # Extract umi sequence
    umi = umi_extractor.extract(reads_in[1].sequence)
    if umi is not None:
        umi_string = f'R0:{umi.start};RX:{umi.sequence};'
    else:
        umi_string = ''
    
    dpm_tag = ''
    if 'TCGAGTCTTGGGTGTTT' in reads_in[1].sequence:
        dpm_tag = 'dp:g'
    if  'CTGACGCTAAGTGCTGAT' in reads_in[1].sequence:
        if len(dpm_tag):
            dpm_tag += 'c'
        else:
            dpm_tag = 'dp:c'
    if len(dpm_tag):
        dpm_tag += ';'
        
    return (
        f'Is:{instrument};{prefix}RN:{run_number};Fc:{flow_cell_id};La:{lane};Ti:{tile};CX:{cluster_x_pos};CY:{cluster_y_pos};{dpm_tag}{barcode_string}{umi_string}'.rstrip(';'),
        even, odd, umi, dpm, y
    )


def fastq_reader(handle: pysam.FastqFile) -> Iterator[pysam.FastxRecord]:
    yield from handle

def recode_reads(fastq_paths: list[str], 
                 output_paths: list[str], 
                 dt: FASTQ_MODALITY_DATATYPES, 
                 library: str,
                 umi_extractor: UmiExtractor,
                 multiplexer_alias: MULTIPLEXER_ALIASES = 'sciT',
                 write_ci:bool = False,
                 
                 ):
    
    
    fastq_in_handles = [pysam.FastqFile(path) for path in fastq_paths] #xopen.xopen(path, 'rt',threads=4) for path in fastq_paths]
    fastq_out_handles = [xopen.xopen(path, 'wt', compresslevel=1,threads=4) for path in output_paths]
    prefix = f'MX:{multiplexer_alias};dt:{dt};LY:{library};' # Cache the prefix of each read
    
    stats = {
        'usable_total':0,
        'unusable_total': 0,
        'missing_umi':0,
        'missing_odd':0,
        'missing_even':0
    }
    for reads_in in zip(*fastq_in_handles):
        
        #verbose = False
        new_header, even, odd, umi, dpm, y = extract_barcodes(reads_in, prefix, umi_extractor, write_ci)
        
        bc_ok = even is not None and odd is not None and umi is not None
        if bc_ok:
            stats['usable_total'] += 1
        else:
            stats['unusable_total'] += 1
            if even is None:
                stats['missing_even'] += 1
            if odd is None:
                stats['missing_odd'] += 1
            if umi is None:
                stats['missing_umi'] += 1

        for i, (record, out_handle) in enumerate(zip(reads_in, fastq_out_handles)):
            record.name = new_header
            # @ optional todo: When the UMI is found, trim the read down? 
            
            out_handle.write(str(record) + '\n')

    return stats


def run():
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Recode sciT demultiplexed (BarcodeIdentification_v1) reads to scmo format"
    )

    argparser.add_argument(
        '--r1_in',
        help="Path to input fastq file for R1",
        type=str, required=True
        )
    argparser.add_argument(
        '--r2_in',
        help="Path to input fastq file for R2",
        type=str
        )
    argparser.add_argument(
        '--r1_out',
        help="Path to output fastq file for R1",
        type=str, required=True
        )
    argparser.add_argument(
        '--r2_out',
        help="Path to output fastq file for R2",
        type=str
        )    
    argparser.add_argument('--dt', help='Datatype', choices=FASTQ_MODALITY_DATATYPES.__args__)
    argparser.add_argument('--multiplexer', help='Multiplexer', choices=MULTIPLEXER_ALIASES.__args__)
    argparser.add_argument('--library', help='Library')
    
    argparser.add_argument('--seq_before_umi', help='Sequence before the UMI', required=True)
    argparser.add_argument('--seq_after_umi', help='Sequence after the UMI', required=True) 
    argparser.add_argument('--expected_umi_length', help='UMIs with lengths different than this value are rejected, when no value is set the length is autodetected', type=int)
    argparser.add_argument('--max_flank_edits', help='Maximum number of edits allowed in each of the two sequences flanking the UMI', type=int, default=1)
    
    argparser.add_argument('--write-ci', help='Write ci field in the BAM file, indicating multiome data', action='store_true')
    argparser.add_argument('--stats', help='Path to yml to write statistics to', type=str)
    argparser.add_argument('-v', '--version', action='version', version=__version__)
    args = argparser.parse_args()
    assert args.r1_in != args.r2_in, "Input files for R1 and R2 must be different"
    assert args.r1_out != args.r2_out, "Output files for R1 and R2 must be different"
    assert args.r2_in is None or args.r2_out is not None, "If R2 is specified, R2_out must be specified"
    
    print(f'sciT_read_recorder version: {__version__}')
    if args.expected_umi_length is not None:
        umi_length = args.expected_umi_length
        print(f"Using expected UMI length {umi_length}")
    else:
        umi_length, confidence = autodetect_umi_len(args.r2_in, args.seq_before_umi, args.seq_after_umi)
        print(f"Using autodetected UMI length {umi_length}, reads with this UMI length is at least {confidence*100:.2f}%")
        
    # Create the UMI extractor
    umi_extractor = UmiExtractor(
        sequence_before_umi=args.seq_before_umi,
        sequence_after_umi=args.seq_after_umi,
        umi_length=umi_length, 
        max_edits=args.max_flank_edits)

    stats = recode_reads(
        fastq_paths = [args.r1_in, args.r2_in],
        output_paths=[args.r1_out, args.r2_out],
        dt=args.dt,
        library=args.library,
        multiplexer_alias = args.multiplexer,
        umi_extractor=umi_extractor,
        write_ci = args.write_ci
    )
    if args.stats is not None:
        with open(args.stats, 'w') as h:
            yaml.dump(stats, h)
    
if __name__ == '__main__':
    run()
    
    