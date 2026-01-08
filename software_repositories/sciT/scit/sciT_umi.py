import random
import re
from math import ceil
from itertools import product
from typing import Iterator
import pysam
import matplotlib.pyplot as plt
import pandas as pd
import re
import numpy as np

def split_sequence_in_parts(sequence, max_edits):
    # Splits the sequence into parts of length max_edits+1, the last part can be shorter when the sequence is not divisible by max_edits+1
    stepper = ceil(len(sequence) / (max_edits+1))
    chunks = list(range(0, len(sequence), stepper))
    parts = [sequence[start: ((start+stepper) if i<len(chunks)-1 else None) ] for i,start in enumerate(chunks)]
    return parts

def get_search_patterns(chunks_hit):
    chunks_miss = [ '.'*len(chunk) for chunk in chunks_hit ]
    search_patterns = []
    for chunk_index, chunk in enumerate(chunks_hit):
        search_patterns.append(  ''.join([ chunk if idx == chunk_index else chunks_miss[idx] for idx in range(len(chunks_hit)) ] ))
    return search_patterns


class SequenceHit:
    def __init__(self, total_sequence: str, hit_sequence: str, start: int, end: int, hamming_dist: int, hamming_dist_before: int, hamming_dist_after: int, matcher):
        self.total_sequence = total_sequence
        self.sequence = hit_sequence
        self.start = start
        self.end = end
        self.total_sequence = total_sequence
        self.hamming_dist = hamming_dist
        self.hamming_dist_before = hamming_dist_before
        self.hamming_dist_after = hamming_dist_after
        
        self.matcher = matcher
        self.umi = self.sequence #[len(self.matcher.sequence_before_umi):len(self.matcher.sequence_before_umi)+self.matcher.umi_length]
        
    def __repr__(self) -> str:
        return f"SequenceHit({self.total_sequence}, umi={self.umi}, hit_sequence={self.sequence}, start={self.start}, end={self.end}, hd={self.hamming_dist})"


class UmiExtractor:
    def __init__(self, sequence_before_umi: str, sequence_after_umi: str, umi_length: int, max_edits: int, MIN_POS: int = 50):
        self.max_edits = max_edits
        self.sequence_before_umi = sequence_before_umi
        self.sequence_after_umi = sequence_after_umi
        patterns_before = get_search_patterns( split_sequence_in_parts(sequence_before_umi, max_edits) )
        patterns_after = get_search_patterns( split_sequence_in_parts(sequence_after_umi, max_edits) )
        self.exact_pattern = sequence_before_umi + '.'*umi_length + sequence_after_umi
        search_patterns = [ pattern_before + '.'*umi_length + pattern_after for pattern_before, pattern_after in product(patterns_before, patterns_after) ]
        self.regex = re.compile( f"(?={'|'.join(search_patterns)})" )
        self.umi_length = umi_length
        self.min_pos = MIN_POS
        
    def find_matches(self, sequence: str) -> Iterator[SequenceHit]:
        for hit in self.regex.finditer(sequence[self.min_pos:]):
            # Check if the hits has less than max edits:
            hit_sequence = sequence[ hit.start() + self.min_pos : hit.start() + self.min_pos + self.umi_length + len(self.sequence_before_umi) + len(self.sequence_after_umi) ]
            
            hamming_dist_before_umi = sum([ ref_char != char for  (char, ref_char) in zip(hit_sequence, self.sequence_before_umi) if ref_char != '.' ])
            if hamming_dist_before_umi > self.max_edits:
                continue
            hamming_dist_after_umi = sum([ ref_char != char for  (char, ref_char) in zip(hit_sequence[self.umi_length+len(self.sequence_before_umi):], self.sequence_after_umi) if ref_char != '.' ])
            if hamming_dist_after_umi > self.max_edits:
                continue

            #hamming_dist = sum([ ref_char != char for  (char, ref_char) in zip(hit.group(), self.exact_pattern) if ref_char != '.' ])
            hamming_dist = hamming_dist_before_umi + hamming_dist_after_umi
            umi = hit_sequence[len(self.sequence_before_umi) : len(self.sequence_before_umi) + self.umi_length]
            yield SequenceHit( sequence, umi, hit.start() + self.min_pos, hit.start()+ self.umi_length + self.min_pos, hamming_dist, 
                                hamming_dist_before_umi, hamming_dist_after_umi,
                                self )
    def extract(self, sequence) -> SequenceHit | None:
        matches = list(self.find_matches(sequence))
        # Select the best match (lowest hamming distance)
        if len(matches) == 0:
            return None
        best_match = min(matches, key=lambda sequence_hit: sequence_hit.hamming_dist)
        assert len(best_match.umi) == self.umi_length, f"UMI length is {len(best_match.umi)} but should be {self.umi_length}"
        return best_match
    

def mutate_sequence(string: str, n_edits: int):
    chars = list(string)
    applied_edits = 0
    while applied_edits < n_edits:
        pos = random.randint(0, len(chars)-1)
        if chars[pos] != 'N':
            chars[pos] = 'N'
            applied_edits += 1
    return ''.join(chars)


def autodetect_umi_len(read_path:str, sequence_before_umi:str, sequence_after_umi:str, sample_n:int = 5_000, MAX_POS = 150, MIN_POS = 50, create_plots=False) -> int | None:

    start_location_obs = pd.Series(np.zeros(MAX_POS))
    umi_len_obs = pd.Series(np.zeros(MAX_POS))

    # First determine the (most common) start location of the UMI
    n_sampled = 0

    pattern_before = re.compile(sequence_before_umi)
    
    skip_due_to_not_demultiplexed = 0
    skip_due_to_MIN_POS = 0
    skip_due_to_MAX_POS = 0
    skip_due_to_no_match = 0

    with pysam.FastqFile(read_path) as h:
        for i,record in enumerate(h):
            #if i%1000==0 and i>0:
            #    print(f"Read {i} reads, skipped, accepted {n_sampled}({100*n_sampled/i:.2f}%), skipped {skip_due_to_not_demultiplexed}({100*skip_due_to_not_demultiplexed/i:.2f}%) due to not demultiplexed, {skip_due_to_MIN_POS} due to MIN_POS, {skip_due_to_MAX_POS} due to MAX_POS, {skip_due_to_no_match}({100*skip_due_to_no_match/i:.2f}%) due to no match", end='\r')
            if record.name.count('[NOT_FOUND]') > 1:
                skip_due_to_not_demultiplexed+=1
                continue
            used_read = False
            for hit in pattern_before.finditer(record.sequence):
                before_umi = hit.start()
                
                if  before_umi < MIN_POS: # or before_umi == -1 
                    skip_due_to_MIN_POS+=1
                    continue
                start = before_umi+len(sequence_before_umi)
                if start>=MAX_POS:
                    skip_due_to_MAX_POS+=1
                    continue
                used_read = True
                start_location_obs[start] += 1

            if used_read:
                n_sampled += 1
                if n_sampled >= sample_n:
                    break
            else:
                skip_due_to_no_match+=1

    # Determine umi length
    start_loc = start_location_obs.idxmax()
    if create_plots:
        start_location_obs.plot()    
        plt.axvline(start_loc, color='red', lw=0.5)
        plt.title(f'Umi usually starts at {start_loc}bp')
        plt.show()

    n_sampled = 0
    before_umi_start = start_loc - len(sequence_before_umi)
    with pysam.FastqFile(read_path) as h:
        for record in h:
            if record.name.count('[NOT_FOUND]') >= 2:
                continue
            
            # Check if the start sequence matches exactly:
            if record.sequence[before_umi_start:start_loc] != sequence_before_umi:
                continue
            
            after_umi = record.sequence[start_loc:].find(sequence_after_umi)
            if after_umi == -1 or after_umi >= MAX_POS:
                continue
            
            umi_len_obs[after_umi] += 1
            n_sampled+=1
            if n_sampled >= sample_n:
                break
    
    umi_len = umi_len_obs.idxmax()
    confidence = umi_len_obs.loc[umi_len] / umi_len_obs.sum()
    
    if create_plots:
        umi_len_obs.plot()
        umi_len = umi_len_obs.idxmax()
        plt.axvline(umi_len, color='red', lw=0.5)
        plt.title(f'Umi length = {umi_len}bp ')
    return umi_len, confidence
