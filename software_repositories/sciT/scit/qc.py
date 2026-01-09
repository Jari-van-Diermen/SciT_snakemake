#!/usr/bin/env python3
import argparse
import pandas as pd
import anndata as ad
import seaborn as  sns
import matplotlib.pyplot as plt
import numpy as np
from singlecellmultiomics.utils import createRowColorDataFrame
import matplotlib.colors
import yaml
import time


def load_loom(loom_path):
    adata = None
    for i in range(1000):
        try:
            adata = ad.read_loom(loom_path)
            break
        except BlockingIOError:
            time.sleep(0.1)
    if adata is None:
        raise Exception(f"Could not open loom file {loom_path}")
    return adata

def loom2cond(loom, condition_keys):
    adata = load_loom(loom)
            
        
    condition_frame = adata.obs[ [k for k in condition_keys if k in adata.obs.columns] ]
    condition_frame.index = loom2df_tx(loom).index
    condition_frame[condition_frame == 'nan'] = np.nan
    return condition_frame
    
def loom2df_tx(loom_path):
    adata = load_loom(loom_path)
    
    df = pd.DataFrame( adata.X.toarray(), 
                        index=prune_multiome_identifier(adata.obs['sample']), 
                        columns=pd.MultiIndex.from_frame(adata.var[['identifier',]])  )
    return df

def prune_multiome_identifier(l):
    return [x.replace('_multiome_','_').replace('_DamID_','_') for x in l]
    
def not_enough_samples_plot():
    
    fig, ax = plt.subplots()
    plt.text(0.5,0.5,'Not enough samples to cluster',horizontalalignment='center')
    sns.despine(left=True, bottom=True)
            
def quality_control(transcriptome_loom_files,
                    transcriptome_clustering_output_path=None,
                    transcriptome_count_threshold = 0,
                    transcriptome_genes_threshold = 1000,
                    config_path = None
                   ):

    sample_qc_status = {}    
    
    if config_path is not None:
        
        try:
            with open(config_path) as h:
                config = yaml.safe_load(h)['conditions']
            condition_keys = config['attributes'] # The attributes to use for condition coloring
        except KeyError:
            print(f'No conditions found in config file {config_path}')
            row_colors=None
            config = None
        
        if config is not None:
            # Parse the colors
            if 'colors' in config:
                attribute_colors = {attribute:{value.replace('\\',''):matplotlib.colors.hex2color(f'#{c}') for value,c in value_colors.items()} for attribute, value_colors in config['colors'].items()}
            else:
                attribute_colors = None
                
            condition_table = pd.concat([loom2cond(p, condition_keys) for p in transcriptome_loom_files]) 
            row_colors, lut = createRowColorDataFrame(condition_table, predeterminedColorMapping=attribute_colors)
            print(row_colors.shape)
            if row_colors.shape[1]==0:
                print('No remaining color information for row colors. No row colors will be used.')
                row_colors = None
    else:
        row_colors=None

    # Transcriptome:
    tx_counts = pd.concat( map(loom2df_tx, transcriptome_loom_files)).fillna(0)
    print(tx_counts)
    sample_qc_status['RNA'] = {}
    for cell, row in tx_counts.iterrows():
        n_reads = row.sum()
        n_genes = (row>0).sum()
        sample_qc_status['RNA'][cell] = (n_reads>=transcriptome_count_threshold) and (n_genes>=transcriptome_genes_threshold)

    tx_counts = tx_counts[pd.Series(sample_qc_status['RNA'])]
    
    cors= pd.DataFrame( 
        np.corrcoef( tx_counts), 
        index=prune_multiome_identifier(tx_counts.index), 
        columns=prune_multiome_identifier(tx_counts.index)
    )
    if cors.shape[0]<=2:
        not_enough_samples_plot()
    else:
        cm = sns.clustermap(
            cors,cmap='RdBu_r',#row_colors=row_colors, col_colors=row_colors, 
            figsize=(11,11),
            dendrogram_ratio=0.05,
            cbar_pos=(-0.05,0.9,0.025,0.1),
            #yticklabels=True,xticklabels=True,
            method='average',
            annot=False,
            row_colors=row_colors,
            col_colors=row_colors
        )
        plt.suptitle(f'Correlation between samples \n Transcriptome',y=1.02)
    if transcriptome_clustering_output_path is not None:
        plt.savefig(transcriptome_clustering_output_path, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
    
    #sns.clustermap(densities,z_score=0,col_cluster=False,row_cluster=False,vmax=3, vmin=-1, figsize=(10,8), cmap='Greys')
    return sample_qc_status

    
def run():
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run quality control on supplied loom files, using provided thresholds"
    )

    ip = argparser.add_argument_group('Inputs')
    ip.add_argument(
        '--tx',
        help="Path to transcriptome loom file, multiple files can be supplied, the selected CONDITION COLUMNS need to be encoded in this file",
        action='append', 
        )
    
    flt = argparser.add_argument_group('Thresholds')
    flt.add_argument('--tx_count_threshold',type=int, default=0, help='Minimum required number of unique RNA reads to pass QC')
    flt.add_argument('--tx_genes_threshold',type=int, default=1000, help='Minimum required number of unique RNA genes to pass QC')
    
    cnd = argparser.add_argument_group('Conditions')
    cnd.add_argument('--config_file',type=str, help="""Path to config.yaml file which has a key conditions:
  attributes: # Attributes used to determine the condition (label) of each sample
  - fusion_construct
  - phase
  colors: # Predetermined colors for attribute values:
    fusion_construct:
      53BP1: 0000ff
    phase: # used fucci reporter colors, colors are HEX without the preceding #
      G1: ff0000
      G1_S: ffcc00
      G2: 00a400""")
    
    out = argparser.add_argument_group('Outputs')

    out.add_argument(
        '--tx_clustering_output_path',
        help="path to transcriptome clustering plot output",
        default=None,
        type=str
        )
    out.add_argument(
        '--qc_csv',
        help="path to qc csv file",
        type=str
        )
    
    args = argparser.parse_args()

    assert args.tx is not None, 'Provide at least one loom file'
    
    sample_qc_status = quality_control(
                    transcriptome_loom_files=args.tx,
                    transcriptome_clustering_output_path=args.tx_clustering_output_path,
                    transcriptome_count_threshold = args.tx_count_threshold,
                    transcriptome_genes_threshold=args.tx_genes_threshold,
                    config_path=args.config_file
    )
    if args.qc_csv is not None:
        qcstat = pd.DataFrame(sample_qc_status)
        qcstat.to_csv(args.qc_csv)
    
if __name__ == '__main__':
    run()
    
    