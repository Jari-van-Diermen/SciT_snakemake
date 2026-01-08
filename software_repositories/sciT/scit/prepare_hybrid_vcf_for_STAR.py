import os
import pysam
import argparse
import pysam
from collections import Counter 

__version__ =  '1'

def prepare_hybrid_vcf(
    vcf_in_path,
    vcf_out_path,
    sample_a,
    sample_b,
    drop_multi_allelic_variants: bool=True # drop variants with more than 2 alleles, eg genotypes like 1|2
):
    """
    Converted a multi sample vcf to a single sample vcf file with phased genotypes
    
    """
    hybrid_sample = f'{sample_a}_{sample_b}'
    rejection_reasons = Counter()

    with pysam.VariantFile(vcf_in_path) as vcf_in:
        header = pysam.VariantHeader()
        # Add a contig (chromosome 1) to our VCF header
        for idx,x in list(vcf_in.header.contigs.items()):
            header.contigs.add(x.name, x.length, name=x.name)
        # Add GT to FORMAT in our VCF header
        header.add_meta('FORMAT', items=[('ID',"GT"), ('Number',1), ('Type','String'),
            ('Description','Genotype')])
        
        header.add_sample(hybrid_sample)
        os.makedirs(os.path.dirname(vcf_out_path), exist_ok=True)
        with pysam.VariantFile(vcf_out_path, header = header, mode = 'w') as vcf_out:
            for record in vcf_in:
                
                alleles_a = set(record.samples[sample_a].alleles)
                alleles_b = set(record.samples[sample_b].alleles)
                if len(alleles_a) != 1  or len(alleles_b) != 1:
                    rejection_reasons['More than one allele']+= 1
                    continue
                a = list(alleles_a)[0]
                b = list(alleles_b)[0]
                if a is None or b is None:
                    rejection_reasons['A or B is none'] += 1
                    continue
                if a == b:# Skip alleles which do not differ from one other
                    rejection_reasons['a=b'] += 1 
                    continue
                # create new record
                clone_basis_record = record
                # recalculate alleles, remove alleles which do not exist in samples
                ref_allele = record.alleles[0]
                alt_alleles = [x for x in record.alleles[1:] if x in [a, b]]
                alleles = tuple([ref_allele] + alt_alleles)
                if drop_multi_allelic_variants and len(alleles) > 2:
                    rejection_reasons['multi-allelic'] += 1
                    continue
                rejection_reasons['pass'] += 1
                # Populate new record
                record_args = {
                    'contig':clone_basis_record.chrom,
                    'start':clone_basis_record.start,
                    'alleles':alleles,
                    'id': clone_basis_record.id,
                    #'info': clone_basis_record.info
                }
                try:
                    new_record = vcf_out.new_record(**record_args)
                except Exception as e:
                    print(record_args)
                    raise e
                new_record.samples[hybrid_sample].alleles = (a, b)
                new_record.samples[hybrid_sample].phased = True
                vcf_out.write(new_record)
    return rejection_reasons

def parse_args():
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Generate hybrid bam file from multi sample vcf file"
    )

    argparser.add_argument(
        '--vcf_in_path',
        help="Multi sample vcf file",
        type=str, required=True
        )
    
    argparser.add_argument(
        '--vcf_out_path',
        help="Hybrid vcf file output path",
        type=str, required=True
        )
    
    argparser.add_argument(
        '--sample_a',
        help="Sample to associate to first haplotype (A)",
        type=str, required=True
        )
  
    argparser.add_argument(
        '--sample_b',
        help="Sample to associate to second haplotype (B)",
        type=str, required=True
        )
    
    argparser.add_argument('-v', '--version', action='version', version=__version__)
    
    return argparser.parse_args()

def run():
    args = parse_args()
    for reason, n in prepare_hybrid_vcf(
        args.vcf_in_path, 
        args.vcf_out_path, 
        args.sample_a, 
        args.sample_b
    ).most_common():
        print(f'{reason}: {n}')
    
if __name__ == '__main__':
    run()