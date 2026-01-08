""" MultiQC module to parse multi-omic outputs """


from typing import Dict
import yaml
import logging
from collections import defaultdict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, heatmap, scatter, table
from multiqc.utils import mqc_colour

import numpy as np
import os

# Initialise the logger
log = logging.getLogger(__name__)
tab20 = ((0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
 (0.6823529411764706, 0.7803921568627451, 0.9098039215686274),
 (1.0, 0.4980392156862745, 0.054901960784313725),
 (1.0, 0.7333333333333333, 0.47058823529411764),
 (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
 (0.596078431372549, 0.8745098039215686, 0.5411764705882353),
 (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
 (1.0, 0.596078431372549, 0.5882352941176471),
 (0.5803921568627451, 0.403921568627451, 0.7411764705882353),
 (0.7725490196078432, 0.6901960784313725, 0.8352941176470589),
 (0.5490196078431373, 0.33725490196078434, 0.29411764705882354),
 (0.7686274509803922, 0.611764705882353, 0.5803921568627451),
 (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
 (0.9686274509803922, 0.7137254901960784, 0.8235294117647058),
 (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
 (0.7803921568627451, 0.7803921568627451, 0.7803921568627451),
 (0.7372549019607844, 0.7411764705882353, 0.13333333333333333),
 (0.8588235294117647, 0.8588235294117647, 0.5529411764705883),
 (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
 (0.6196078431372549, 0.8549019607843137, 0.8980392156862745))

def format_cell_name(cell_name):
    return cell_name.replace('_multiome','').replace('_DamID_','_')

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Single Cell Multi Omics",
            anchor="scmo",
            href="barbansonbiotech.com",
            info="SCMO module for DAMID"
        )
        
         # Load sample whitelists if available
        self.show_cells = set()
        n_sample_files = 0
        for f in self.find_log_files('scmo/samplelist'):
            n =0 
            n_sample_files+=1
            for cell in f['f'].splitlines():
                cell = cell.replace('_multiome','').replace('_DamID_','_') # For brevity this multiome string is left out in the QC report
                self.show_cells.add(cell)
                n+=1
            log.info(f"Loaded {n} cell names from {f['fn']}")
        log.info(' '.join(list(self.show_cells)[:10]))
        
                
        self.general_header = dict()
        self.meta_data = dict()
        self.meta_header = dict()
        # parse meta data file
        # meta information can be library or sample/cell specific.
        # for bulk protocols with little samples the meta can be produced for each sample
        # for single cell samples this is obviously too much data
        for f in self.find_log_files("scmo/meta"):
            meta = self.parse_meta(f)
            self.meta_data.update( meta )
            del meta
        #log.info(self.meta_data.keys())
        #log.info(self.meta_data['SLEB-ES-DamLBR-i3_1'])
        self.damid_count_threshold = 5000
        self.transcriptome_count_threshold = 5000
        self.transcriptome_genes_threshold = 5000
        self.show_allele_specific_data = None # Will be autodetected
        self.cell_min_damid_count = 0
        self.cell_min_transcriptome_count = 0
        
        self.report_top_samples = 100
        self.report_bottom_samples = 100
       

        if n_sample_files == 0:
            self.show_cells = None
            
        if os.path.exists('config.yaml'):
            with open('config.yaml') as h:
                
                config = yaml.safe_load(h)
                if 'DamID' in config and 'scidam' in config['DamID']:
                    self.cell_min_damid_count = config['DamID']['scidam'].get('cell_min_damid_count', self.cell_min_damid_count)
                    self.cell_min_transcriptome_count = config['DamID']['scidam'].get('cell_min_transcriptome_count', self.cell_min_transcriptome_count)

                self.attribute_colors = {}
                
                qc_config = config.get('QC', {})

                self.damid_count_threshold = qc_config.get('damid_count_threshold', self.damid_count_threshold)
                self.transcriptome_genes_threshold = qc_config.get('transcriptome_genes_threshold', self.transcriptome_genes_threshold)
                self.transcriptome_count_threshold = qc_config.get('transcriptome_count_threshold', self.transcriptome_count_threshold)
                
                self.report_top_samples = qc_config.get('report_top_samples', self.report_top_samples)
                self.report_bottom_samples = qc_config.get('report_bottom_samples', self.report_bottom_samples)
                
                if 'conditions' in config:
                    config= config['conditions']
                    self.condition_keys = config['attributes'] # The attributes to use for condition coloring
                    if 'colors' in config:
                        self.attribute_colors = {
                            attribute:{value.replace('\\',''):
                                f'#{c}'
                                for value,c in value_colors.items()} 
                            for attribute, value_colors in config['colors'].items()}

                log.info(self.attribute_colors)
                        


        # Load tagging information
        """
        'overall_statistics':overall_statistics,
        'datatype_statistics':datatype_statistics,
        'per_data_type':{
            'DamID':dam_statistics,
            'RNA':transcriptome_statistics
        }}
        
        G9c-64nM_1:
            duplicate: 1
            qcfail: 25416
            unique: 6
        
        """
        
        self.tagging_data = dict()
        self.tagging_header = dict()
        self.removed_cells_tagging_data = {}
        max_show_in_log = 10
        shown_in_log = 0
        for f in self.find_log_files('scmo/tagging'):
            yml = yaml.safe_load( f["f"] )
            library = f['s_name'].replace(' | se_tagging','').replace(' | tagging','')
            log.info(f'Parsing tag file {f["fn"]}')
            
            for datatype, dt_data in yml['per_data_type'].items():
                for sample, sample_info in dt_data.items():
                    sample = format_cell_name(sample)
                    if sample.startswith('Unknown_'): #Skip reads with no assigned sample
                        continue

                    if shown_in_log<max_show_in_log:
                        log.info(f'Keeping data from sample {sample} {library} {datatype}, is ignored: {self.show_cells is not None and sample not in self.show_cells}')
                        shown_in_log+=1
                    if self.show_cells is not None and sample not in self.show_cells:
                        self.add_sample_info_to_tagging_data(target_dict=self.removed_cells_tagging_data, sample=sample,  library=library, datatype=datatype, sample_info=sample_info)
                        continue

                    self.add_sample_info_to_tagging_data(target_dict=self.tagging_data, sample_info=sample_info, sample=sample,  library=library, datatype=datatype)
                    if shown_in_log<max_show_in_log:
                        log.info(f'Keeping data from sample {sample} {library} {datatype}')
                        shown_in_log+=1
        
        scidamdemux_data = dict()
        for f in self.find_log_files('scmo/scidamdemux'):
            # Load to yaml:
            yml = yaml.safe_load( f["f"] )
            library = f['s_name'].split(' | ')[0] # f['s_name'] Will be something like "MB-Dam-10K | scidam.demux.RNA" or "MB-Dam-10K | scidam.demux.DamID"
            total = yml['unusable_total'] + yml['usable_total']
            if not library in scidamdemux_data:
                scidamdemux_data[library] = dict()
            if f['s_name'].endswith('.RNA'):
                scidamdemux_data[library]['RNA_total_input_reads'] = total
                scidamdemux_data[library]['RNA_total_demultiplexed_reads'] = yml['usable_total']
                log.info( "SAMPLE NAME FOR RNA FILE:" + f['s_name'] ) 
            else:
                scidamdemux_data[library]['DamID_total_input_reads'] = total
                scidamdemux_data[library]['DamID_total_demultiplexed_reads'] = yml['usable_total']
                log.info( "SAMPLE NAME FOR DAMID FILE:" + f['s_name'] ) 
                
            log.info( "SAMPLE NAME FOR TAG FILE:" + library ) #.replace(' | se_tagging','').replace(' | tagging','')
            

        
        # Remove ignored samples:
        self.tagging_data = self.ignore_samples(self.tagging_data)
       
        # Add the statistics to the meta table
        # calculate library aggregates too:
        # Add allele specific statistics too
        """
                
        `mapped_allele_`
        `qcfail_allele_{read_haplotype}`
        `unique_allele_{read_haplotype}`
        `duplicate_allele_{read_haplotype}`

        """
        self.allele_specific_keys = [
            'DamID_mapped_allele_1',
            'DamID_mapped_allele_2',
            'DamID_unique_allele_1',
            'DamID_unique_allele_2',
            'DamID_qcfail_allele_1',
            'DamID_qcfail_allele_2',
            'RNA_mapped_allele_1',
            'RNA_mapped_allele_2',
            'RNA_unique_allele_1',
            'RNA_unique_allele_2',
            'RNA_unique_genes_allele_1',
            'RNA_unique_genes_allele_2',
        ]
        
        
        found_allele_specific_data = False
        library_aggregates = dict()
        for library, library_data in (self.removed_cells_tagging_data | self.tagging_data).items():
            if not library in library_aggregates:
                library_aggregates[library] = { 
                                               'DamID_reads':0,
                                               'DamID_unique':0, 
                                               'DamID_QC':0,
                                               'RNA_reads':0, 
                                               'RNA_unique':0, 
                                               'RNA_unique_genes':0,
                                               'RNA_QC':0,
                                               'All_modalities_QC':0
                                               }
                # Add allele specific keys
                for k in self.allele_specific_keys:
                    library_aggregates[library][k] = 0
                
            log.info(f'Adding data for {library}')
            
            for sample,sample_data in library_data.items():
                if '_multiome_' in sample:
                    sample = sample.replace('_multiome','')
                elif '_DamID_' in sample:
                    sample = sample.replace('_DamID','')
                
                write_meta = True
                if self.show_cells is None: # DamID logic
                    if not sample in self.meta_data:
                        continue
                else:
                    
                    if not sample in self.meta_data:
                        if sample in self.show_cells: # Only add the sample if we want to show it in the statistics table
                            #log.info(f'No meta data provided for {sample}, dropping this sample/cell') # For normal DamID we skipped these cells

                            self.meta_data[sample] = {}
                        else:
                            write_meta =False
                        
                if write_meta:
                    self.meta_data[sample]['Library'] = library
                
                for key, value in sample_data.items():
                    #log.info(f'Adding {key} to {sample}')
                    if not found_allele_specific_data and key in self.allele_specific_keys:
                        found_allele_specific_data = True
                    
                    if write_meta:
                        self.meta_data[sample][key] = value
                    if key in library_aggregates[library]:
                        library_aggregates[library][key]+=value
                        
                all_pass =True
                some_pass = False
                if sample_data.get('DamID_unique',0)>=self.damid_count_threshold:
                    library_aggregates[library]['DamID_QC'] += 1
                    some_pass = True
                else:
                    all_pass=False
                if (sample_data.get('RNA_unique',0)>=self.transcriptome_count_threshold) and (sample_data.get('RNA_unique_genes',0)>=self.transcriptome_genes_threshold):
                    library_aggregates[library]['RNA_QC'] += 1
                    some_pass = True
                else:
                    all_pass=False
                    
                if  all_pass:
                    library_aggregates[library]['All_modalities_QC'] +=1
                if some_pass and not write_meta:
                    log.warning(f'Cell {sample} passes QC but will not be written to the table')
                # 
        
        for lib, d in scidamdemux_data.items():
            if not lib in library_aggregates:
                library_aggregates[lib] = {}
            library_aggregates[lib].update(d)
            
        
        if self.show_allele_specific_data is None and found_allele_specific_data:
            self.show_allele_specific_data = True
        ## Calculate percentages
        for sample, sample_data in self.meta_data.items():
            
            for name, a,b in (('DamID Mapping %', 'DamID_mapped','DamID_reads'),
                              ('DamID Usable %', 'DamID_unique','DamID_reads'),
                              ('RNA Mapping %', 'RNA_mapped','RNA_reads'),
                              ('RNA Usable %', 'RNA_unique','RNA_reads'),
                              
                              ):
                if a not in sample_data or b not in sample_data:
                    continue
                self.meta_data[sample][name] = (sample_data[a]/sample_data[b])*100
        
        ## Calculate percentages for library aggregates:
        for sample, sample_data in library_aggregates.items():
            
            for name, a,b in (('DamID Mapping %', 'DamID_mapped','DamID_reads'),
                              ('DamID Usable %', 'DamID_unique','DamID_reads'),
                              ('RNA_Mapping_pct', 'RNA_mapped','RNA_reads'),
                              ('RNA Usable %', 'RNA_unique','RNA_reads'),
                              ('DamID Demux %', 'RNA_total_demultiplexed_reads','RNA_total_input_reads'),
                              ('RNA Demux %', 'RNA_total_demultiplexed_reads','RNA_total_input_reads'),
                              ):
                if a not in sample_data or b not in sample_data or sample_data[b] == 0:
                    continue
                library_aggregates[sample][name] = (sample_data[a]/sample_data[b])*100

        self.update_meta_headers(self.meta_data)
        
        # Write the collected data
        self.write_data_file(self.meta_data, "scmo-sample-information-table")
        
        #### Library table
        t_config = {
            "id": "scmo-libraries",
            "namespace": "library_information",
            "title": "Library_information",
            "no_violin": True,
            "save_file": True
            
        }
        
        self.general_stats_addcols(library_aggregates, self.general_header)

        log.info("Adding library table")
        log.info(library_aggregates)
        log.info("Header")
        log.info(self.general_header)
        p = table.plot(data=library_aggregates, headers=self.general_header, pconfig=t_config)
        log.info(p)

        self.add_section(
            name="Library statistics",
            anchor="lib-stats",
            description="Statistics aggregated per library",
            plot=p
            )
        #####
        
        ### SC table:
        self.add_meta_table()
        
        
        
        # Create plate overview data
        cols384 = [ f'{x}' for x in range(1,25)] 
        rows384 = [x for  x in 'ABCDEFGHIJKLMNOP']
        colormaps = {
            
            'kw' : [[0, "#000000"], [1, "#ffffff"]],
            'wr' : [[0, "#000000"], [1, "#ff0000"]],
            
        }
        for metric, pct, description, minval, maxval, colormapname in ( 
                                     ('DamID_unique', None, 'Number of reads with a DamID adapter, correct cut site motif and not duplicate',0, 10000,'kw' ),
                                     ('RNA_unique', None, 'Number of reads with a RNA adapter, assigned to a gene and not duplicate',0, 10000,'kw' ),
                                     ('DamID_qcfail', 'DamID_reads', 'Percentage of reads with a DamID adapter (vs total reads with DamID adapter), but not assigned to a cut site',0, 100,'wr' ),
                                     ('RNA_qcfail', 'RNA_reads', 'Percentage of reads with a RNA adapter (vs total reads with RNA adapter), not assigned to a gene',0, 100,'wr' )):
            

            htmls = []
            for library, library_data in self.tagging_data.items():
                log.info(f'Generating plate plot for {library} {metric}, {len(library_data)} cells')
                plate = np.zeros( (16, 24) )
                for sample,sample_data in library_data.items(): 
                    # Check if the well is available:
                    sample = sample.replace('_multiome','') #.replace('_DamID','')
                    
                    if sample not in self.meta_data or not 'Well' in self.meta_data[sample]:
                        #log.info(f'{sample} missing in meta_data')
                        continue
                    #else:
                    #log.info(f'{sample} found in meta_data')
                    
                    col = self.meta_data[sample]['Well'][1:]
                    row = self.meta_data[sample]['Well'][0]
                    # log.info(self.meta_data[sample]['well'], row, col)
                    if not row in rows384:
                        log.info(f'{row} not in {rows384}')
                        row, col = col, row
                    
                    if metric in sample_data:
                        value = sample_data[metric]
                        
                        if pct is not None and not pct in sample_data:
                            log.info(f'Tried to calculate {metric} over {pct}, but a value for {pct} is missing for sample {sample}, library: {library}')
                            continue
                        if pct is not None:
                            normalizer = sample_data[pct]
                            value = (value/normalizer)*100
                        plate[rows384.index(row),cols384.index(col)] = value

                if plate.sum().sum()>0:
                    log.info(f'PLATE VIEW GENERATED {library}-{metric}')
                    htmls.append(heatmap.plot(plate.tolist(), 
                                             xcats=[f'{col} ' for col in cols384],  # The space behind col works around a bug where the first two x-ticks overlap
                                             ycats=rows384,
                                             pconfig={
                                                'title':f'{library} {metric}',
                                                'display_values':False,
                                                'xlab':'column',
                                                'ylab':'row',
                                                'zlab':metric,
                                                'min':minval,
                                                'max':maxval,
                                                "xcats_samples": False, 
                                                "ycats_samples": False,
                                                "colstops": colormaps[colormapname],
                                                'id':f'plate-view-{library}-{metric}'
                                             }
                                             ))

                else:
                    log.info(f'No data available for {library} {metric}')
            if len(htmls)>0:
                log.info(htmls[0])
                log.info(str(htmls[0]))
                self.add_section(
                    name=f"Plate view {metric}",
                    anchor=f"plate-view-{metric}",
                    description=description,
                    content = ''.join([x.interactive_plot() for x in htmls])
                )

        sdat = self.create_scatterplot_data(self.tagging_data, x_attr='DamID_unique', y_attr='RNA_unique', force_label=None)
        
            

        if len(self.removed_cells_tagging_data)>0:
            
            for label, data in self.create_scatterplot_data(self.removed_cells_tagging_data, 
                                                              x_attr='DamID_unique', 
                                                              y_attr='RNA_unique', 
                                                 fixed_color='rgba(120,120,120,0.2)', 
                                                 force_label='removed'
                                                 ).items():
                
                if label not in sdat:
                    log.info(f"Adding {label} to scatterplot data, there are no threshold passing cells")
                    sdat[label] = []
                log.info(f"Writing {len(data)} (removed) data points for {label}, already added: {len(sdat[label])}")
                sdat[label] += data
        else:
            log.warn("No removed cells (self.removed_cells_tagging_data), removed labels will not be present")

        for label, data in sdat.items():
            log.info(f"Ended up with {len(data)} data points for {label}")
            log.info(f'Eg: {data[:10]}')
            
        self.add_section(
            name='Transcriptome vs Dam reads per sample',
            anchor="dam_vs_rna",
            description="Transcriptome vs Dam reads per sample",
            plot=scatter.plot(sdat,
                              { 
                               #"square": True,
                                "logswitch": True, 
                                "xlog": True,           
                                "ylog": True,
                                "xlab":'DamID unique reads',
                                "ylab":'RNA unique reads',
                                "id":'transcriptome_vs_damid',
                                'title':"Transcriptome vs Dam reads per sample",
                                # 'ymax':1_000_000_000,
                                # 'ymin':1_000_000_000,
                                # 'xmin':1,
                                # 'xmax':1,
                              }
                              )
        )
        

    def add_sample_info_to_tagging_data(self, target_dict, sample_info:Dict[str, int | str], sample:str , library: str, datatype:str ): 
    
        for key, value in sample_info.items():
            target_key = f'{datatype}_{key}' # Go from a nested to a flat data structure so for example RNA:{'reads':0} to {RNA_reads:0}
            if not library in target_dict:
                target_dict[library] = dict()
            if not sample in target_dict[library]:
                target_dict[library][sample] = dict()
                
            #log.info(f"{library} dt:{datatype}, sample:{sample}, {target_key}:{value}")
            if target_key in target_dict[library][sample]:
                # Try to add the numbers 
                
                try:
                    target_dict[library][sample][target_key] += value
                except Exception as e:
                    ##log.info(f"{e}, dt:{datatype}, sample:{sample}, {target_key}:{value}")
                    target_dict[library][sample][target_key] = value
            else:
                target_dict[library][sample][target_key] = value
                
                    
    def create_scatterplot_data(self, data, x_attr, y_attr, fixed_color=None, force_label=None):
        sdata = {}
        min_value = 1
        for i,(library, library_data) in enumerate(data.items()):
            if force_label is None:
                label = library
            else:
                label = force_label
            log.info(f"Creating scatterplot data for {library}, will get label {label}")
            if label not in sdata:
                sdata[label] = []
            for sample, sample_data in library_data.items(): # Sample is usually a cell here
                x,y = sample_data.get(x_attr,0.1), sample_data.get(y_attr,0.1)
                x = x if x is not None and x>min_value else min_value
                y = y if y is not None and y>min_value else min_value
                sdata[label].append({
                    'x':x,
                    'y':y,
                    'color': ( 
                        ('rgba(' + ','.join([str(int(c*255))  for c in tab20[i%20]] ) + ',0.8)') 
                        if fixed_color is None else fixed_color ),
                })
        return sdata
        
    def add_meta_table(self):
        t_config = {
            "id": "scmo-sample-information-table",
            "namespace": "scmo",
            "title": "Sample meta",
            "no_violin": True,
            
        }
        # Reduce number of cells in this table by using self.report_bottom_samples and self.report_top_samples thresholds
        sort_by = 'DamID_unique'
        top_samples = [sample for sample, d in sorted( [ (sample, d.get(sort_by,0)) for sample, d in self.meta_data.items()], key= lambda x: x[1], reverse=True)[:self.report_top_samples]]
        bottom_samples = [sample for sample, d in sorted( [ (sample, d.get(sort_by,0)) for sample, d in self.meta_data.items()], key= lambda x: x[1], reverse=False)[:self.report_bottom_samples]]
        select_samples = set(top_samples + bottom_samples)
        
        description="Sample statistics and meta information"
        if len(select_samples) < len(self.meta_data):
            description += f" only the top {self.report_top_samples} and bottom {self.report_bottom_samples} samples shown, from a total of {len(self.meta_data)} samples" 
            
        
        self.add_section(
            name="Sample information",
            anchor="scmo-sample-information-table-section",
            description=description,
            plot=table.plot({sample:d for sample, d in self.meta_data.items() if sample in select_samples },  self.meta_header, t_config),
        )


    def update_meta_headers(self, meta):
        default_hidden_columns = set([
            'DamID_mapped',
            'DamID_qcfail',
            'DamID_unmapped',
            'DamID_duplicate',
            'RNA_mapped',
            'RNA_qcfail',
            'RNA_unmapped',
            'RNA_duplicate'
            
        ])

                                             

        general_cols = set(['DamID_unique', 'DamID_reads', 'RNA_unique', 'RNA_reads',
                            'DamID_QC','RNA_QC','All_modalities_QC',
                            'DamID_total_input_reads',
                            'DamID_total_demultiplexed_reads',
                            'RNA_total_input_reads',
                            'RNA_total_demultiplexed_reads',
                            'DamID Demux %', 'DamID Mapping %', 'DamID Usable %',  
                            'RNA Demux %',  'RNA_Mapping %', 'RNA Usable %',
                            
                            
                            ])

        descriptions = {
            
            'RNA_total_input_reads':'SciDAM only: Number of reads which are sent to the demultiplexer (after trimming)',
            'RNA_total_demultiplexed_reads': 'SciDAM only: Number of reads assigned to a sample by the demultiplexer',
            'DamID_total_input_reads':'SciDAM only: Number of reads which are sent to the demultiplexer (after trimming)',
            'DamID_total_demultiplexed_reads': 'SciDAM only: Number of reads assigned to a sample by the demultiplexer',
            'DamID_reads':'Number of read pairs assigned to this sample and the DamID modality by the demultiplexer, note that for scidam this includes reads where the sample barcode or umi is missing.',
            'DamID_mapped':'Number of read pairs assigned to this sample and the DamID modality by the demultiplexer, and mapped to the reference genome',
            'RNA_mapped':'Number of read pairs assigned to this sample and the transcriptome by the demultiplexer, and mapped to the reference genome',
            'RNA_reads':'Number of reads assigned to this sample and the transcriptome by the demultiplexer',
            'RNA_unique_genes':f'Number unique genes found in the sample QC pass (green) when more than {self.transcriptome_genes_threshold} unique genes are detected ',
            'DamID_unique':f'Number of reads which are correctly located on a cut site and have a unique molecular identifier and cell barcode. QC pass (green) when more than {self.damid_count_threshold} reads',
            'DamID_qcfail':'Number of reads which cannot be assigned to a cut site, a unique molecular identifier and cell barcode need to be present.',
            'RNA_unique':f'Number of reads which can be assigned to a gene and have a unique molecular identifier and cell barcode. QC pass (green) when more than {self.transcriptome_count_threshold} reads ',
            'RNA_qcfail':'Number of reads which cannot be assigned to a gene',
            
            'DamID Mapping %':'Percentage of multiplexed reads which map to the reference genome',
            'RNA Mapping %':'Percentage of multiplexed reads which map to the reference genome',
            'DamID Usable %':'Percentage of multiplexed reads which map to a cut site and have a unique UMI/Cell barcode',
            'RNA Usable %':'Percentage of multiplexed reads which map to the reference genome, can be assigned to a gene and have a unique UMI/Cell barcode',
            'RNA Demux %':'Percentage of reads which can be assigned a barcode and umi (scidam only)',
            'DamID Demux %':'Percentage of reads which can be assigned a barcode and umi (scidam only)'
            
        }
        
        if self.show_allele_specific_data:
            descriptions.update({
                'DamID_mapped_allele_1':'Number of read pairs assigned to this sample and the DamID modality by the demultiplexer, mapped to the reference genome, and uniquely assignable to allele 1.',
                'DamID_mapped_allele_2':'Number of read pairs assigned to this sample and the DamID modality by the demultiplexer, mapped to the reference genome, and uniquely assignable to allele 2.',
                'DamID_unique_allele_1':'Number of allele 1 reads which are correctly located on a cut site and have a unique molecular identifier and cell barcode.',
                'DamID_unique_allele_2':'Number of allele 2 reads which are correctly located on a cut site and have a unique molecular identifier and cell barcode.',
                'DamID_qcfail_allele_1':'Number of allele 1 reads which cannot be assigned to a cut site',
                'DamID_qcfail_allele_2':'Number of allele 2 reads which cannot be assigned to a cut site',
                'RNA_mapped_allele_1': 'Number of RNA reads which are assigned to allele 1',
                'RNA_mapped_allele_2': 'Number of RNA reads which are assigned to allele 2',
                'RNA_unique_allele_1': 'Number of unique RNA reads which are assigned to allele 1',
                'RNA_unique_allele_2': 'Number of unique RNA reads which are assigned to allele 2',
                'RNA_unique_genes_allele_1': 'Number of unique genes assigned to allele 1',
                'RNA_unique_genes_allele_2': 'Number of unique genes assigned to allele 2'
            })
            
            
            
        # These override the default descriptions (The thresholds for QC pass are not applicable)
        general_descriptions = {
            'DamID_unique':f'Number of reads which are correctly located on a cut site and have a unique molecular identifier and cell barcode.',
            'RNA_unique':f'Number of reads which can be assigned to a gene and have a unique molecular identifier and cell barcode.',
            'DamID_QC':'Number of cells/samples which pass the DamID quality control',
            'RNA_QC':'Number of cells/samples which pass the transcriptome quality control',
            'All_modalities_QC': 'Number of cells/samples which pass the DamID and transcriptome quality control',
            'RNA_total_input_reads': "Number of reads in the library before demultiplexing",
            'RNA_total_demultiplexed_reads': "Number of reads with usable barcode",
            'DamID_total_input_reads': "Number of reads in the library before demultiplexing",
            'DamID_total_demultiplexed_reads': "Number of reads with usable barcode"
        }
        if self.show_allele_specific_data:
            general_descriptions.update({
                'DamID_unique_allele_1':'Number of unique DamID reads assocated to allele 1',
                'DamID_unique_allele_2':'Number of unique DamID reads assocated to allele 2',
                'RNA_unique_allele_1':'Number of unique RNA reads assocated to allele 1',
                'RNA_unique_allele_2':'Number of unique RNA reads assocated to allele 2',
                'DamID_mapped_allele_1':'Number of mapped DamID reads assocated to allele 1',
                'DamID_mapped_allele_2':'Number of mapped DamID reads assocated to allele 2',
                'RNA_mapped_allele_1':'Number of mapped RNA reads assocated to allele 1',
                'RNA_mapped_allele_2':'Number of mapped RNA reads assocated to allele 2'
                })

        integer_columns = ['DamID_unique','RNA_unique','DamID_QC','RNA_QC','All_modalities_QC','DamID_reads','DamID_mapped','RNA_reads','RNA_unique_genes','DamID_unique','DamID_qcfail','RNA_unique','RNA_qcfail',
                           'RNA_total_input_reads','RNA_total_demultiplexed_reads', 'DamID_total_demultiplexed_reads', 'DamID_total_input_reads', 'DamID_total_input_reads'
            
        ]
        
        for sample, sample_data in meta.items():
            
            for key,values in sample_data.items():
                if not key in self.meta_header:
                    
                    header_data = {
                        "title": key.replace('_',' '),
                        "description": descriptions.get(key,key),
                        "hidden": key in default_hidden_columns
                    }
                    if key in integer_columns:
                        header_data[ "format" ] = "{:,.0f}" 
                    
                    non_colored_data = header_data.copy()
                    
                    # Perform coloring by attribute
                    if key in self.attribute_colors:
                        rules = dict()
                        for value, color in self.attribute_colors[key].items():
                            rules[value] = [{'s_eq':value},]
                            
                        u =  {"cond_formatting_rules":rules,
                              "cond_formatting_colours":[
                             {k:str(v)} for k,v in self.attribute_colors[key].items()]
                         }
                        log.info(u)
                        header_data.update(u)
                    
  
                    if key in ('DamID_unique','RNA_unique'):
                        header_data.update({
                            "scale": "RdYlGn",
                            "min":0,
                            "max":10000,
                           
                        })
                    if key == 'DamID_unique':
                        header_data.update({
                            "cond_formatting_rules": {
                            "pass": [ {"gt": str(self.damid_count_threshold)},],   
                            "fail": [{"lt": str(self.damid_count_threshold)}]
                           }}
                        )    
                    if key == 'RNA_unique': 
                        header_data.update({
                            "cond_formatting_rules": {
                            "pass": [ {"gt": str(self.transcriptome_count_threshold-1)},],   
                            "fail": [{"lt": str(self.transcriptome_count_threshold)}]
                           }}
                        )    

                    if key == 'RNA_unique_genes': 
                        header_data.update({

                            "cond_formatting_rules": {
                            "pass": [ {"gt": str(self.transcriptome_genes_threshold-1)},],   
                            "fail": [{"lt": str(self.transcriptome_genes_threshold)}]
                           }}
                        )  

                    self.meta_header[key] = header_data
                    #### add to general columns:
                    if key in general_cols:
                        self.general_header[key] = non_colored_data | {"hidden":False, 
                                                                       'title':key.replace('_', ' '),
                                                                       "description":general_descriptions.get(key, descriptions.get(key,key)),
                                                                       }
                    ####
        for key, desc in general_descriptions.items():
            if not key in self.general_header:
                self.general_header[key] = {'title':key.replace('_', ' '), 'description':desc}
                if key in integer_columns:
                    self.general_header[key][ "format" ] = "{:,.0f}" 
                    

        
    def parse_meta(self, f):
        yml = yaml.safe_load( f["f"] )
        library = f['s_name'].replace(' | meta','')
        library_meta = dict()
        sample_meta = dict()
        # Load library information:
        ignore_keys = set(['single-end-fastq-files',
                       'paired-end-fastq-files',
                       'cell-meta',
                       'multiomes'])
        for key,val in yml.items():
            if key in ignore_keys:
                continue
            key = sanitize_key(key)
            library_meta[key] = val
        
        if 'library-meta' in yml and self.show_cells is not None:
            
            # Add all the meta to all known samples
            
            cells_for_lib = [cell for cell in self.show_cells if cell.rsplit('_',1)[0] == library]
            
            for sample in cells_for_lib:
                if not sample in sample_meta:
                    sample_meta[sample] = dict()
                
            
            for key, value in yml['library-meta'].items():
                key = sanitize_key(key)
                library_meta[key] = value
                
                for sample in cells_for_lib:
                    sample_meta[sample][key] = value
                

        if 'cell-meta' in yml and 'meta' in yml['cell-meta']:
            #assert 'primary-index' in yml['cell-meta'], f'{library} does not have a primary index in cell-meta'
            #primary_index = yml['cell-meta']['primary-index']
            
            for indexer, cell_data in yml['cell-meta']['meta'].items():
                # Sanitize the dict
                # Remap all caps variation of wEll to Well
                cell_data = sanitize_dict({k:v for k,v in cell_data.items()})
                                           
                index = f'{library}_{indexer}'
                sample_meta[index]=dict()
                sample_meta[index].update(library_meta)
                sample_meta[index].update(cell_data)
        return sample_meta

def sanitize_key(k):
    if type(k) == str:
        if k.lower()=='well':
            return 'Well'
        return k.replace('#','_n')
    else:
        return k

def sanitize_dict(d):
    return {sanitize_key(k):v for k,v in d.items()}
        