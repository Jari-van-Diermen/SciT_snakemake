#!/usr/bin/env python3
import argparse
import loompy
import pandas as pd
import pysam 
from scit.sciT_tagger import parse_ds_field
from singlecellmultiomics.bamProcessing import get_samples_from_bam
import numpy as np
from singlecellmultiomics.features import FeatureContainer
from multiprocessing import Pool
from singlecellmultiomics.utils import pool_wrapper
from scit import wasp_read_is_assigned_to_haplotype_1, wasp_read_is_assigned_to_haplotype_2

__version__ = '2'

def cell_name_to_attributes(cell_name:str, sample_to_meta:dict) -> dict:
    # Good idea: load meta data from meta files into the loom file 
    try:
        library, cell = cell_name.split('_',1)
    except ValueError:
        print('Could not obtain cell and library (missing underscore)')
        print(cell_name)
        raise
    return {
        'library':library,
        'cell_index':cell,
        'sample':cell_name,
        'CellID':cell_name
    } | sample_to_meta.get(cell_name, dict())
    
def iter_gene_data(featurecontainer):
    # Just converts the data field of the featurecontainer to dict and yields it
    for contig,start,stop,gene_id,strand,data in featurecontainer:
        yield contig,start,stop,gene_id,strand, dict(data)

def count_features(bam_path: str, contig: str, feature_tag: str, sample_order: list[str]) -> dict[str, np.ndarray]:
    feature_counts = {}
    with pysam.AlignmentFile(bam_path,threads=3) as h:
        for read in h.fetch(contig):
            if read.is_qcfail or read.is_duplicate or read.is_secondary or read.is_supplementary:
                continue
            feature = read.get_tag(feature_tag)
            if feature not in feature_counts:
                feature_counts[feature] = np.zeros(len(sample_order))
            sample = read.get_tag('SM')
            feature_counts[feature][sample_order[sample]] += 1
    return feature_counts

def assign_read_to(read: pysam.AlignedSegment) -> tuple[str, str]:
    if wasp_read_is_assigned_to_haplotype_1(read):
        return 'allele_a', ''
    if wasp_read_is_assigned_to_haplotype_2(read):
        return 'allele_b', ''
    return 'allele_ambiguous', ''

def empty_feature_counts(allele_specific = True) -> dict[str, dict[str, np.ndarray]]:
    if allele_specific:
        return {'allele_a':{}, 'allele_b':{}, 'allele_ambiguous':{}, '':{}}
    return {'':{}}

def generic_single_end_tx_read_pass_function(read: pysam.AlignedSegment) -> bool:
    if read.is_qcfail or read.is_duplicate or read.is_secondary or read.is_supplementary or not read.is_mapped:
        return False
    return True

def count_features_allelic(bam_path, contig, feature_tag, sample_order) -> dict[str, dict[str, np.ndarray]]:
    feature_counts = empty_feature_counts()
    with pysam.AlignmentFile(bam_path,threads=3) as h:
        for read in h.fetch(contig):
            if not generic_single_end_tx_read_pass_function(read):
                continue
            feature = read.get_tag(feature_tag)
            for assign_to in assign_read_to(read):
                if feature not in feature_counts[assign_to]:
                    feature_counts[assign_to][feature] = np.zeros(len(sample_order))
                sample = read.get_tag('SM')
                feature_counts[assign_to][feature][sample_order[sample]] += 1
            
    return feature_counts


def run():
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Extract an loom file from a tagged BAM file"
    )

    argparser.add_argument(
        'bam_in',
        help="Path to input bam files",
        type=str
        )
    argparser.add_argument(
        '--loom_out',
        help="Write loom file here",
        type=str
        )
    
    argparser.add_argument(
        '--gene_tag',
        help="Tag which determines the gene",
        type=str,
        default='GN'
        )
    
    argparser.add_argument(
        '--annotations',
        help="Annotation file, when supplied this is used to also encode the GENE NAME in the loom file ",
        type=str
        )
    
    argparser.add_argument(
        '--allele_specific',
        help="Generate allele specific counts",
        action='store_true'
        )
    
    argparser.add_argument(
        '--threads',
        help="Worker threads",
        type=int
        )
    argparser.add_argument('-v', '--version', action='version', version=__version__)
        
    args = argparser.parse_args()
    feature_counts = empty_feature_counts(args.allele_specific)
    
    feature_tag = args.gene_tag
    with pysam.AlignmentFile(args.bam_in,threads=3) as h:
        sample_order = {sample:i for i,sample in enumerate(list( get_samples_from_bam(h) ))}
        sample_to_meta = dict()
        for rg in h.header['RG']:
            meta = parse_ds_field(rg['DS'])
            sample_to_meta[ rg['SM'] ] = meta
        
        # Obtain sample order
        cmds = [
            ((count_features_allelic if args.allele_specific else count_features ),{
                
                'bam_path':args.bam_in, 
                'contig':c, 
                'feature_tag':feature_tag,
                'sample_order':sample_order
            })
            for c in h.references
        ]
        with Pool(args.threads) as workers:
            for i,r in enumerate(workers.imap_unordered(pool_wrapper, cmds)):
                if args.allele_specific:
                    for layer, counts in r.items():
                        feature_counts[layer].update(counts)
                else: 
                    feature_counts[''].update(r)
                print(f'{len(cmds)-i} contigs remaining ', end='\r')
        print('Counting finished!    ')

    
    # Add genes without any reads
    if args.annotations is not None:
        print('Loading annotations')
        gn = FeatureContainer()
        gn.loadGTF(args.annotations,store_all=True)
        
        for layer in feature_counts:
            missing = set( (gene_id for contig,start,stop,gene_id,strand,data in gn )).difference( set(feature_counts[layer].keys()))
            if len(missing)>0:
                print(f'({layer}) Adding {len(missing)} genes with zero counts')
                print(f'Some examples: {",".join(list(missing)[:10])}')
                for m in missing:
                    feature_counts[layer][m] = np.zeros(len(sample_order))
                
    print('Creating loom file')
    # Convert dictionaries to dataframes and sort by sample names
    for layer in list(feature_counts.keys()): 
        feature_counts[layer] = pd.DataFrame(feature_counts[layer], 
                                             index=sample_order.keys()).sort_index().sort_index(axis=1).T
        print(f'Layer "{layer}" size is : {feature_counts[layer].shape}')

    # Drop samples without valid reads
    #feature_counts = feature_counts.loc[ :, feature_counts.sum(1)>0 ]
    if args.annotations is None:
        print('No annotations provided, just using the gene id as names without more information')
        row_attributes = pd.DataFrame({col:{'Gene':col, 'Name':col} for col in feature_counts[''].index}).T.to_dict("list")
    else:
        
        
        id_to_atts = {
            gene_id:{'Gene':gene_id, 
                     'identifier':data.get('gene_name', gene_id),
                     'name':data.get('gene_name', None), 
                     'contig':contig, 
                     'start':start, 
                     'stop':stop, 
                     'strand':strand, 
                     'gene_version':data.get('gene_version'), 
                     'gene_biotype':data.get('gene_biotype')}
            for contig,start,stop,gene_id,strand,data in iter_gene_data(gn)}
        
        row_attributes = pd.DataFrame({col:{'Gene':col} | id_to_atts.get(col, {})  for col in feature_counts[''].index}).T.to_dict("list")
        
    col_attributes = pd.DataFrame( {str(cell):cell_name_to_attributes(cell, sample_to_meta)  for cell in feature_counts['']}).T.to_dict("list")
    loompy.create(args.loom_out, 
              {name: layer.values for name, layer in feature_counts.items()},
              row_attributes,
              col_attributes
             )
    
    
if __name__ == '__main__':
    run()
    
    