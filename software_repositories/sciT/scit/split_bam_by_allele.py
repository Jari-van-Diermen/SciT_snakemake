#!/usr/bin/env python3
import argparse
import pysam


def split_alignments_by_tag(bam_in, tag, value_to_target_path, missing_tag_out_path=None):
    with pysam.AlignmentFile(bam_in, 'rb', threads=2) as bam_in_file:
        output_handles = {
            value: pysam.AlignmentFile(path, header=bam_in_file.header, mode='wb')
            for value, path in value_to_target_path.items()
        }
        if missing_tag_out_path is not None:
            missing_handle = pysam.AlignmentFile(missing_tag_out_path, header=bam_in_file.header, mode='wb')

        for record in bam_in_file:
            if not record.has_tag(tag):
                if missing_tag_out_path is not None:
                    missing_handle.write(record)
                continue
            tag_value = record.get_tag(tag)
            if tag_value not in output_handles:
                continue
            output_handles[tag_value].write(record)
        for handle in output_handles.values():
            handle.close()
        for path in value_to_target_path.values():
            pysam.index(path)
        if missing_tag_out_path is not None:
            missing_handle.close()
            pysam.index(missing_tag_out_path)
            

def run():
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Split bam file by HP tag"
    )

    argparser.add_argument(
        '--bam_in',
        help="Path to input bam file",
        type=str,required=True
        )
    argparser.add_argument(
        '--hp1_out',
        help="Path to output bam containing allele/haplotype 1 reads",
        type=str,required=True
        )
    argparser.add_argument(
        '--hp2_out',
        help="Path to output bam containing allele/haplotype 2 reads",
        type=str,required=True
        )
    argparser.add_argument(
        '--ambiguous_out',
        help="Path to output bam containing ambigous reads",
        type=str
        )
    args = argparser.parse_args()
    split_alignments_by_tag(bam_in=args.bam_in, tag='HP', missing_tag_out_path=args.ambiguous_out , value_to_target_path={
        1: args.hp1_out,
        2: args.hp2_out,
    })