#!/usr/bin/env python3
import argparse
import datetime
from multiprocessing import Pool
import sys
from typing import Any, Dict, Hashable, List, Literal, Set, Tuple
from singlecellmultiomics.universalBamTagger import QueryNameFlagger
import pysam
from contextlib import ExitStack
from singlecellmultiomics.bamProcessing import merge_bams, sorted_bam_file, write_program_tag
from collections import defaultdict, Counter
from singlecellmultiomics.utils import pool_wrapper
import yaml
from yaml.representer import Representer
import yaml
from uuid import uuid4
from scit import wasp_read_is_assigned_to_haplotype_1, wasp_read_is_assigned_to_haplotype_2
import matplotlib.pyplot as plt
import logging
from typing import Iterator

__version__ = '15'

READ_GROUP_REJECTION_REASON = Literal['BC',] | None # Missing data caused the read to be rejected (not assigned to a sample)
LogLevelType = Literal[logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL]
log_level_dict = {
    'DEBUG': logging.DEBUG,
    'INFO': logging.INFO,
    'WARNING': logging.WARNING,
    'ERROR': logging.ERROR,
    'CRITICAL': logging.CRITICAL
}

def iter_pairs(path: str, contig: str | None = None, threads:int=3) -> Iterator[list[pysam.AlignedSegment]]: 
    """Iterate over pairs of reads in an !> UNSORTED <! BAM file. If contig is specified, only iterates over reads in that contig."""
    records = []
    read_name = None
    if contig is not None:
        with pysam.AlignmentFile(path, threads=threads) as h:
            for i,record in enumerate(h.fetch(contig)):
                if record.is_secondary or record.is_supplementary:
                    continue
                    
                if read_name is not None and record.query_name!=read_name:
                    records.sort(key=lambda x: not x.is_read1)
                    yield records
                    records = []
                records.append(record)
                read_name = record.query_name 
            if len(records):
                # Sort to read1 first:
                yield records
    else:
        with pysam.AlignmentFile(path, threads=threads) as h:
            for i,record in enumerate(h):
                if record.is_secondary or record.is_supplementary:
                    continue
                    
                if read_name is not None and record.query_name!=read_name:
                    records.sort(key=lambda x: not x.is_read1)
                    yield records
                    records = []
                records.append(record)
                read_name = record.query_name 
            if len(records):
                yield records

def recurse_add(source: Dict[Hashable, Any], target: Dict[Hashable, Any]):
    # Recursively add the nested source dictionary to the target
    for key,value in source.items():
        if isinstance(value, dict):
            if key not in target:
                target[key] = Counter()
                #print('Creating',key)
            recurse_add(value, target[key])
        else:
            #print(key, value)
            try:
                target[key] += value
            except KeyError:
                target[key] = value

def prepare_multiome_translation_table(multiomes_path: str) -> Tuple[
                                                            Dict[Tuple[str, ...], Any], # origin_key (tags) -> (tag, value)
                                                            List[str]
                                                            ]: 
    with open(multiomes_path,'r') as h:
        d = yaml.safe_load(h)
    # First obtain all tags used in the origin
    origin_tags: set[str] = set()
    for mapper in d:
        origin_tags.update(set(mapper['origin'].keys()))
    # Fix the order of these tags, this order will determine the key order in the translation table
    origin_tags_list = list(origin_tags)
    translation_table: Dict[Tuple[str, ...], Any] = dict()
    for mapper in d:
        origin_key = tuple( str(mapper['origin'].get(tag)) for tag in origin_tags_list)
        translation_table[origin_key] = mapper['target']
    return translation_table, origin_tags_list

    
def apply_read_group(reads: List[pysam.AlignedSegment], 
                     multiome_indices_match: bool = False) -> Tuple[str, Dict[str,str], READ_GROUP_REJECTION_REASON]: 
    
    lib: str = reads[0].get_tag('LY')
    dt: str = reads[0].get_tag('dt')
    
    # When the ci tag is set this means the origin cell is known.
    # If not it means the data is singleton
    rejection_reason: READ_GROUP_REJECTION_REASON = None
    if not reads[0].has_tag('BC'): # This can happen for scidam data
        bc: str = 'NotFound' 
        rejection_reason = 'BC'
    else:
        bc: str = reads[0].get_tag('BC') # for normal DamID data BC is always set
    
    if not reads[0].has_tag('bi'): # For normal DamID data the bi tag is always set
        sample = f"Unknown_{lib}"
        description = f'BC:{bc};dt:{dt};multiome:0'
    else:
        if reads[0].has_tag('ci'): # A cell index is set, this is multiome data
            bi: int = reads[0].get_tag('bi')  # Barcode index
            indexer: int = reads[0].get_tag('ci') 
            sample = f'{lib}_multiome_{indexer}'
            description = f'BC:{bc};ci:{indexer};bi:{bi};dt:{dt};multiome:1'
        elif multiome_indices_match:
            # The cell index matches the barcode
            bi: int = reads[0].get_tag('bi')
            sample = f'{lib}_multiome_{bi}'
            description = f'BC:{bc};ci:{bi};bi:{bi};dt:{dt};multiome:1'
        else:
            indexer = reads[0].get_tag('bi') 
            sample = f'{lib}_{dt}_{indexer}'
            description = f'BC:{bc};bi:{indexer};dt:{dt};multiome:0'
        
    pu = dt + '.'+ reads[0].get_tag('Fc') + '.' + reads[0].get_tag('La') + f'.{bc}.' + sample 
    id = pu
    for read in reads:
        read.set_tag('RG', id)
        read.set_tag('SM', sample)
        
    return id, {'ID':id, 
                'LB':lib, 
                'PL':'Illumina',
                'SM':sample, 
                'PU':pu, 
                'BC':bc,
                'DS':description
                }, rejection_reason


                
def _parse_single_ds_kv_pair(kvpair: Tuple[str, str]) -> Tuple[str, str]:
    return tuple( (kv.replace(':','-colonreplace-').replace(';','-semicolonreplace-')  for kv in kvpair) ) # type: ignore | Iterates first over the key, then the value
    
def parse_ds_field(s: str) -> Dict[str, str]:
    # Parse a ds field  made of a key-value string "k:v;k2:v2" to a dict
    parts = s.split(';' )
    return dict( _parse_single_ds_kv_pair(kv.split(':',1)) for kv in parts ) # type: ignore


def tag_damid_file(file_in_path: str,
                   damid_out_path: str,
                   cut_site_assigner: CutSiteAssigner | None, 
                   transcriptome_out_path: str | None = None,
                   translation_table: Dict[Tuple[str], Any] | None = None, #  , # origin_key (tags) -> (tag, value)                   
                   translation_tags: List[Set[str]] | None = None,
                   meta: dict | None = None,
                   contig: str | None = None,
                   transcriptome_only: bool = False,
                   gene_tag: str = 'GN', # Tag used to assign a read to a gene
                   multiome_indices_match: bool = False,
                   min_mapping_quality: int = 0, 
                   ):

    base_tagger = QueryNameFlagger()
    write_transcriptome = transcriptome_out_path is not None
    damid_molecules_seen = set() # mi,contig,pos
    transcriptome_molecules_seen = defaultdict(set) # sample -> {gene, RX}

    # { 'ID': ID, 'LB':library,
    #             'PL':platform,
    #             'SM':sampleLib,
    #             'PU':readGroup }
    damid_read_groups = dict()
    transcriptome_read_groups = dict()
    
    overall_statistics = defaultdict(Counter) # cell->{total_reads,mapped_reads,unmapped_reads,} etc
    datatype_statistics = defaultdict(Counter)
    dam_statistics = defaultdict(Counter)
    transcriptome_statistics = defaultdict(Counter) 
    with ExitStack() as exitstack:

        infile = exitstack.enter_context(pysam.AlignmentFile(file_in_path))
        
        out_header = infile.header.as_dict()
        
        write_program_tag(
            out_header,
            program_name='damid_tagger',
            command_line= " ".join(
                sys.argv),
            version= __version__,
            description=f'DamID tagger, executed at {datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')


        if not transcriptome_only:
            damid_out = exitstack.enter_context( sorted_bam_file(damid_out_path, header=out_header, read_groups=damid_read_groups) )
        if write_transcriptome:
            transcriptome_out = exitstack.enter_context( sorted_bam_file(transcriptome_out_path, header=out_header,read_groups=transcriptome_read_groups) )
        
        for records in iter_pairs(file_in_path, contig=contig):
            base_tagger.digest(records) # This reads the query name and extracts tags set by the demultiplexer
            # Pick a read which is not None, preferably R1 if available
            tagged_read = pick_tagged_read(records)
            apply_multiome_mapping(translation_table, translation_tags, records, tagged_read) # Apply multiome mapping (connects multiple modalities through ci tag)
            read_haplotype = get_haplotype_information(records, tagged_read) # Apply and obtain haplotype information
            rg_id, rg_def, read_group_rejection_reason = apply_read_group(records, multiome_indices_match) # Try to get a read group, fails when no barcode(s) available
            sample_key = tagged_read.get_tag('SM')
            assert type(sample_key) == str, f"Sample key (SM) is not a string: {sample_key}"
            update_overall_statistics(overall_statistics, tagged_read, sample_key, read_group_rejection_reason)
            # Check the datatype of the read:
            try:
                dt = tagged_read.get_tag('dt')
                if dt is None: # The datatype cannot be determined
                    continue
            except KeyError:
                if transcriptome_only:
                    continue
                dt = 'DamID' # Fall back to DamID
            datatype_statistics[sample_key][dt] +=1 
            
            if dt=='RNA':
                if not write_transcriptome:
                    continue
                if rg_id not in transcriptome_read_groups:
                    transcriptome_read_groups[rg_id] = rg_def
                
                is_duplicate, is_qcfail = qc_and_deduplicate_transcriptomic_read(gene_tag, transcriptome_molecules_seen, transcriptome_out, records, tagged_read, sample_key, read_haplotype)
                update_transcriptome_statistics(transcriptome_statistics, tagged_read, read_haplotype, sample_key, is_duplicate, is_qcfail)
                    
            else:
                cut_site_prediction = cut_site_assigner.predict(*records)
                is_valid_read =  cut_site_prediction.is_valid
                assert cut_site_prediction.cut_site_informative_read is not None, "Cut site prediction has no associated informative read"
                if not cut_site_prediction.cut_site_informative_read.has_tag('BC') or not cut_site_prediction.cut_site_informative_read.has_tag('RX'): # Not enough barcoding information available.
                    is_valid_read = False
                    
                if cut_site_prediction.cut_site_informative_read.is_unmapped or cut_site_prediction.cut_site_informative_read.mapping_quality < min_mapping_quality:
                    is_valid_read = False
                
                if not is_valid_read:
                    is_duplicate = False
                else:
                    molecule_hash = (cut_site_prediction.cut_site_informative_read.reference_name, 
                                     cut_site_prediction.predicted_cut_position, 
                                     cut_site_prediction.cut_site_informative_read.get_tag('MI') 
                                        if cut_site_prediction.cut_site_informative_read.has_tag('MI') 
                                        else cut_site_prediction.cut_site_informative_read.has_tag('BC')+cut_site_prediction.cut_site_informative_read.has_tag('RX'), 
                                     read_haplotype)
                    
                    if molecule_hash in damid_molecules_seen:
                        is_duplicate = True
                    else:
                        is_duplicate = False
                        damid_molecules_seen.add(molecule_hash)
                        
                #STATISTICS UPDATE:
                update_dam_statistics(dam_statistics, cut_site_prediction.cut_site_informative_read, read_haplotype, sample_key, is_duplicate, is_valid_read)
                write_read_duplicate_and_rg(damid_read_groups, damid_out, records, rg_id, rg_def, is_duplicate)
                
        # Write unique genes seen:
        for sample_key, molecules in transcriptome_molecules_seen.items():
            transcriptome_statistics[sample_key]['unique_genes']  += len(set( (gene_id for gene_id, umi, haplotype in molecules)))
            # And haplotype specific:
            for allele in (1,2):
                transcriptome_statistics[sample_key][f'unique_genes_allele_{allele}']  += len(set( (gene_id for gene_id, umi, haplotype in molecules if haplotype==allele) ))
        
        # Add meta information to the read groups if such meta data is available
        if meta is not None:
            read_group_meta = extract_read_group_meta(meta, damid_read_groups)
            for rg_id, attributes in read_group_meta.items():
                # Append the collected meta to the DS field
                damid_read_groups[rg_id]['DS']  += ';' + ( ';'.join([f'{k}:{v}' for k,v in attributes.items()]) )
                
            
    if not transcriptome_only:
        # plx.simple_bar(list(observed_offsets_fwd.keys()), 
        #             list(observed_offsets_fwd.values()), 
        #             width = 100, title = 'Observed offsets after correction fwd ')
        # plx.show()
        
        # plx.simple_bar(list(observed_offsets_rev.keys()), 
        #             list(observed_offsets_rev.values()), 
        #             width = 100, title = 'Observed offsets after correction rev ')
        # plx.show()
        pass
            
    return {
        'overall_statistics':overall_statistics, 
        'datatype_statistics':datatype_statistics, # cell-> dt -> n
        'per_data_type':{
            'DamID':dam_statistics,
            'RNA':transcriptome_statistics
        }}

def pick_tagged_read(records):
    tagged_read = records[0]
    for read in records:
        if read is None:
            continue
        if read.is_read1:
            tagged_read = read
            break
        elif read.is_read2 and tagged_read is None:
            tagged_read = read
    return tagged_read
    # per_data_type:
    # DamID: {}
    # RNA:
    #     MB-Dam-10K_multiome_0:
    #     duplicate: 237
    #     mapped: 1787
    #     qcfail: 2291
    #     reads: 3085
    #     unique: 557
    #     unique_genes: 555
    #     unique_genes_allele_1: 0
    #     unique_genes_allele_2: 0
    #     unmapped: 1298

def write_read_duplicate_and_rg(damid_read_groups, damid_out, records, rg_id, rg_def, is_duplicate):
    for read in records:
        if is_duplicate:
            read.is_duplicate = True
        damid_out.write(read)
        if rg_id not in damid_read_groups:
            damid_read_groups[rg_id] = rg_def

def qc_and_deduplicate_transcriptomic_read(gene_tag:str, 
                                           transcriptome_molecules_seen:Set, 
                                           transcriptome_out:pysam.AlignmentFile, 
                                           records:List[pysam.AlignedSegment], 
                                           tagged_read: pysam.AlignedSegment, 
                                           sample_key:str, 
                                           read_haplotype: Literal[1,2] | None):
    is_duplicate = False
    is_qcfail = False
    if not tagged_read.has_tag('RX') or not tagged_read.has_tag(gene_tag) or tagged_read.get_tag(gene_tag)=='-' or tagged_read.get_tag(gene_tag).startswith('__') or not tagged_read.has_tag('BC'):
                    # this is a qcfail molecule
        is_qcfail = True                        
    else:
        gene_id = tagged_read.get_tag(gene_tag)
        molecule_hash = (gene_id, tagged_read.get_tag('RX'), read_haplotype)
                    
        if molecule_hash in transcriptome_molecules_seen[sample_key]:
            is_duplicate = True
        else:
            is_duplicate = False
            transcriptome_molecules_seen[sample_key].add(molecule_hash)
                    
    for tagged_read in records:
        if is_duplicate:
            tagged_read.is_duplicate = True
        if is_qcfail:
            tagged_read.is_qcfail = True 
        transcriptome_out.write(tagged_read)
    return is_duplicate,is_qcfail


def update_overall_statistics(overall_statistics, tagged_read, sample_key, read_group_rejection_reason):
    overall_statistics[sample_key]['reads'] +=1 
    if tagged_read.is_mapped:
        overall_statistics[sample_key]['r1_mapped'] +=1 # @TODO: the counting is not necessary based on R1 anymore (just the name is incorrect)
    else:
        overall_statistics[sample_key]['r1_unmapped'] +=1 # @TODO: the counting is not necessary based on R1 anymore
    if read_group_rejection_reason is not None:
        overall_statistics[sample_key]['not_demultiplexed'] +=1

def apply_multiome_mapping(translation_table, translation_tags, records, read):
    origin_key = tuple( (str(read.get_tag(tag)) if read.has_tag(tag) else None) for tag in translation_tags )
    if origin_key in translation_table:
                # Apply:
        for tag,value in translation_table[origin_key].items():
            for record in records:
                record.set_tag(tag,value)

def get_haplotype_information(records: List[pysam.AlignedSegment], tagged_read: pysam.AlignedSegment) -> Literal[1,2] | None:
    if tagged_read.has_tag('vW'):
        if wasp_read_is_assigned_to_haplotype_1(tagged_read):
            read_haplotype = 1
            for r in records:
                r.set_tag('HP',1)
        elif wasp_read_is_assigned_to_haplotype_2(tagged_read):
            read_haplotype = 2
            for r in records:
                r.set_tag('HP',2)
        else:
            read_haplotype = None    
    else:
        read_haplotype = None
    return read_haplotype

def update_transcriptome_statistics(transcriptome_statistics, 
                                    tagged_read: pysam.AlignedSegment, 
                                    read_haplotype: Literal[1,2] | None,
                                    sample_key: str, 
                                    is_duplicate: bool, 
                                    is_qcfail: bool):
    transcriptome_statistics[sample_key]['reads'] += 1
    if tagged_read.is_mapped:
        transcriptome_statistics[sample_key]['mapped'] += 1
    else:
        transcriptome_statistics[sample_key]['unmapped'] += 1
                    
    if is_qcfail:
        transcriptome_statistics[sample_key]['qcfail'] +=1
    elif not is_duplicate:
        transcriptome_statistics[sample_key]['unique'] +=1
    else:
        transcriptome_statistics[sample_key]['duplicate'] +=1
                    
                # Allele specific
    if read_haplotype is not None:
        transcriptome_statistics[sample_key][f'mapped_allele_{read_haplotype}'] += 1
        if is_qcfail:
            transcriptome_statistics[sample_key][f'qcfail_allele_{read_haplotype}'] +=1
        elif not is_duplicate:
            transcriptome_statistics[sample_key][f'unique_allele_{read_haplotype}'] +=1
        else:
            transcriptome_statistics[sample_key][f'duplicate_allele_{read_haplotype}'] +=1

def update_dam_statistics(dam_statistics: dict, 
                          informative_read: pysam.AlignedSegment, 
                          read_haplotype: Literal[1,2] | None, 
                          sample_key: str, 
                          is_duplicate: bool, 
                          is_valid_cut: bool):
    dam_statistics[sample_key]['reads'] += 1
    if informative_read.is_mapped:
        dam_statistics[sample_key]['mapped'] += 1
    else:
        dam_statistics[sample_key]['unmapped'] += 1
    if not is_valid_cut:
        dam_statistics[sample_key]['qcfail'] +=1
    elif not is_duplicate:
        dam_statistics[sample_key]['unique'] +=1
    else:
        dam_statistics[sample_key]['duplicate'] +=1            
    # Allele specific
    if read_haplotype is not None:
        dam_statistics[sample_key][f'mapped_allele_{read_haplotype}'] += 1
        if not is_valid_cut:
            dam_statistics[sample_key][f'qcfail_allele_{read_haplotype}'] +=1
        elif not is_duplicate:
            dam_statistics[sample_key][f'unique_allele_{read_haplotype}'] +=1
        else:
            dam_statistics[sample_key][f'duplicate_allele_{read_haplotype}'] +=1

def extract_read_group_meta(meta, damid_read_groups) -> Dict[str, Dict[str,Any]]: 
    dt_indexer = meta['cell-meta']['primary-index']
    # Create a mapping which connects the indexer (bi tag) to a list of readgroups
    indexer_to_readgroups = defaultdict(list)
    for rgid, read_group_data in damid_read_groups.items():
        rg_info = parse_ds_field( read_group_data.get('DS') )
        if rg_info['dt']!=dt_indexer:
            continue
                # Some read groups will not have a bi flag (scidam only) when no (complete) BC is present
        if not 'bi' in indexer_to_readgroups:
            continue
        indexer_to_readgroups[ rg_info['bi'] ].append( rgid )
    logging.info('Cell identifier to read groups:')
    for indexer, readgroups in indexer_to_readgroups.items():
        logging.info(f' {indexer} -> {",".join(readgroups)}')
        
    read_group_meta = defaultdict(dict) # Read group -> meta to add
    # Obtain library wide meta information and apply it to all read groups, it can be overriden by cell specific meta:
    if 'library-meta' in meta:
        for k,v in meta['library-meta'].items():
            v = str(v).replace(':','-colonreplace-').replace(';','-semicolonreplace-')                
            k = str(k).replace(':','-colonreplace-').replace(';','-semicolonreplace-')
            for rgid, read_group_data in damid_read_groups.items():
                read_group_meta[rgid][k] = v
                
    if 'meta' in meta['cell-meta']:
        for index, kvpairs in  meta['cell-meta']['meta'].items():
            for k,v in kvpairs.items():
                v = str(v).replace(':','-colonreplace-').replace(';','-semicolonreplace-')                
                k = str(k).replace(':','-colonreplace-').replace(';','-semicolonreplace-')                
                for associated_rg in indexer_to_readgroups[str(index)]:
                    read_group_meta[associated_rg][k] = v

    return read_group_meta
    

def run():
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Assign cut sites and deduplicate a DamID library. Make sure the input file is sorted by query name! The outputs are sorted by coordinate"
    )

    argparser.add_argument(
        '-bam_in',
        help="Path to input bam file (Can be mixture of transcriptome and Dam)",
        type=str,required=True
        )
    argparser.add_argument(
        '-damid_out',
        help="Path to output DamID bam file",
        type=str,required=False
        )
    argparser.add_argument(
        '-transcriptome_out',
        help="Path to output transcriptome bam file, when only this flag is set the data is assumed to be transcriptome only",
        type=str,required=False
        )
    argparser.add_argument(
        '-reference_fasta_path',
        help="Path to reference fasta file",
        type=str,required=True, action='append'
        )
    argparser.add_argument(
        '-statistics',
        help="Path to statistics yaml",
        type=str,required=False
        )
    argparser.add_argument(
        '-motif',
        default='GATC',
        type=str,required=False, help='Cut site motif'
        )
    
    argparser.add_argument(
        '-multiomes',
        type=str,required=False, help='Multiomes file (yaml)'
        )
    argparser.add_argument(
        '-meta',
        type=str,required=False, help='Meta file (yaml), should contain a meta key, primary indexer and then key: attribute: value'
        )
    
    argparser.add_argument(
        '--threads',
        default=1,
        type=int, help='Multithread by splitting by chromosome, this can only be done when the input is single end and coordinate sorted'
        )
    
    argparser.add_argument(
        '-gene_tag',
        type=str,default='GN',help='Tag present in the input bam file which associates reads to a gene, for example GN for STAR, and XF for HTSEQ'
        )
    
    argparser.add_argument(
        '-min_mapping_quality',
        help="Minimum mapping quality for a mapped read to be taken into account",
        type=int, default=10
        )
    
        
    argparser.add_argument(
        '--scidam',
        help="The data is scidam, enables scidam based distance to cut-site function",
        action='store_true'
        )
    
    argparser.add_argument(
        '--pdf_threshold',
        help="SciDAM only: The minimum probability for a DamID cut to be considered (this probability is stored in pC tag), set to 0 to disable this filter",
        type=float, default=0.0005
        )
    argparser.add_argument(
        '-cut_density_plot',
        help="Path to density plot",
        type=str,required=False
        )
    
    argparser.add_argument(
        '--log-level',
        choices=list(log_level_dict.keys()),
        default='INFO')
        
    
    argparser.add_argument(
        '-v',
        action='version',
        version=__version__
        )
    
    
    args = argparser.parse_args()
    
    # Set logging based on args:
    logging.basicConfig(
        level=args.log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    logging.info("Starting damid-tagger")
    logging.info(f"Arguments: {args}")
    translation_table, translation_tags = prepare_multiome_translation_table(args.multiomes)
    
    if args.damid_out is None:
        transcriptome_only= True
    else:    
        transcriptome_only= False

    if not transcriptome_only:
        logging.info("Determining cut locations in refererence")
        centered_coord_sorted_cut_sites = CenteredCoordSortedCutSites(fasta_path= args.reference_fasta_path, 
                                                                      motif="GATC"
                                                                    )
        # For scidam use:
        if args.scidam:
            logging.info('sciDAM mode')
            cut_site_assigner = CutSiteAssigner_SciDam_Multi(
               centered_coord_sorted_cut_sites = centered_coord_sorted_cut_sites,
               dist_clip_after= 1000,
               dist_clip_before= 1000,
               p_threshold=args.pdf_threshold
               )
            cut_site_assigner.fit(args.bam_in, sample_n = 1_000_000)
        else:
            logging.info('DamID mode')
            cut_site_assigner = CutSiteAssigner_DamID(
                centered_coord_sorted_cut_sites = centered_coord_sorted_cut_sites
                )
            cut_site_assigner.fit(args.bam_in, sample_n = 25_000)
            # plx.simple_bar(list(cut_distance_distribution_fwd.keys()), 
            #             list(cut_distance_distribution_fwd.values()), 
            #             width = 100, title = 'Motif distance from FWD read start')
            # plx.show()
            
            # plx.simple_bar(list(cut_distance_distribution_rev.keys()), 
            #             list(cut_distance_distribution_rev.values()), 
            #             width = 100, title = 'Motif distance from REV read start')
            # plx.show()
            
            # Obtain max index:

            logging.info(f'Using fixed offset of {cut_site_assigner.offset_fwd} for fwd strand')
            logging.info(f'Using fixed offset of {cut_site_assigner.offset_rev} for rev strand')
    else:
        cut_site_assigner = None # Not needed for transcriptome
    
    
    if not transcriptome_only and args.cut_density_plot is not None and cut_site_assigner is not None:
        # Write plot with probability density of p cut relative to the read end
        cut_site_assigner.plot()
        plt.savefig(args.cut_density_plot)
    
    if args.meta is not None:
        with open(args.meta) as h:
            meta = yaml.safe_load(h)
    else:
        meta = None
    
    multiome_indices_match = False
    if meta is not None:
        if meta['type'] == 'scDamT':# check if data is multiome
            if 'multiomes' not in meta: # Check if indices are provided
                multiome_indices_match = True
                logging.warning("Multiome data found without provided indices, assuming cell barcode indices are matching between transcriptome and genome")
    logging.info("Writing results")    
    
    if args.threads>1:
        stats = defaultdict(Counter)
        with Pool(args.threads) as workers, pysam.AlignmentFile(args.bam_in) as h:
            contigs = list(h.references) + ['*']
            temp_folder = './temp'
            cmds = []
            for contig in contigs:
                temp_bam_path_dam = f'{temp_folder}/{uuid4()}.dam.bam'
                temp_bam_path_tx = f'{temp_folder}/{uuid4()}.tx.bam'
                
                cmds.append((
                    tag_damid_file,
                    {
                        'file_in_path': args.bam_in,
                        'cut_site_assigner': cut_site_assigner,
                        'damid_out_path': temp_bam_path_dam,
                        'transcriptome_out_path': temp_bam_path_tx,
                        'translation_table' : translation_table,
                        'translation_tags' : translation_tags,
                        'meta': meta,
                        'contig' : contig,
                        'transcriptome_only':transcriptome_only,
                        'gene_tag':args.gene_tag,
                        'multiome_indices_match':multiome_indices_match,
                        'min_mapping_quality':args.min_mapping_quality
                    })
                )
            
            # Execute:
            for result in workers.imap(pool_wrapper,cmds):
                # Merge the statistics
                recurse_add(result, stats)
            
            if not transcriptome_only: 
                merge_bams( [cmd['damid_out_path'] for _,cmd in cmds], args.damid_out, threads=args.threads)
            merge_bams( [cmd['transcriptome_out_path'] for _,cmd in cmds], args.transcriptome_out, threads=args.threads)
    else:
        stats = tag_damid_file(args.bam_in,
                    args.damid_out,
                    cut_site_assigner,
                    transcriptome_out_path=args.transcriptome_out,
                    translation_table=translation_table, 
                    translation_tags=translation_tags,
                    meta= meta,
                    transcriptome_only=transcriptome_only,
                    multiome_indices_match=multiome_indices_match,
                    min_mapping_quality=args.min_mapping_quality
                    )
    

    if args.statistics is not None:
        # Interpret defaultdict and Counter as a dictionary
        yaml.add_representer(defaultdict, Representer.represent_dict)
        yaml.add_representer(Counter, Representer.represent_dict)
        if not transcriptome_only:
            pass
            # stats['offsets'] = {
            #     #'distribution_fwd':cut_distance_distribution_fwd,
            #     #'distribution_rev':cut_distance_distribution_rev,
                
            #     # 'fwd':cut_site_assigner.fwd_offset,
            #     # 'rev':cut_site_assigner.rev_offset (Not valid for scidam)
            # }
        with open(args.statistics,'w') as h:
            yaml.dump(stats,h)
    
    logging.info('Done')
    
    
if __name__ == '__main__':
    run()
    
    