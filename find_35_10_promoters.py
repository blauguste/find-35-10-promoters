from Bio import SeqIO, SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import nt_search
import csv

#gb_path = 'NC_004350_2.gb'
#ref_record = SeqIO.read(open(gb_path), 'genbank')
#ref_sequence = ref_record.seq
ref_sequence = Seq('CCCTTGACACCCCCCCCCCCCCCCCCCTATAATCCCCCCCCCCGGGGGGGGGGCCCCCCCCCCGGGGGGGGGGGGGGGG', IUPAC.unambiguous_dna)
print(ref_sequence)
print(len(ref_sequence))

# List of transcription start sites. 
### Csv with [0] uid, [1] threshold, [2] genome, [3] feature name, [4] strand, [5] TSS
TSS_list = []
with open('sample_tss_list.csv', 'r') as TSS_input:
    TSS_reader = csv.reader(TSS_input, delimiter = ',')
    next(TSS_reader, None) # ignore the header row
    for row in TSS_reader:
        TSS_list.append((row[0], row[1], row[2], row[3], row[4], int(row[5])))
print(TSS_list)

# Promoter sequences to look for: 
### In -10 promoters: A2 and T6 always conserved plus 2 out of 4 variable sites (total of 4/6).
### In -35 promoters: T1 and T2 always conserved plus 2 out of 4 variable sites (total of 4/6).
ten_promoters = ['TATNNT', 'TANANT', 'TANNAT', 'NATANT', 'NATNAT', 'NANAAT']
thirtyfive_promoters = ['TTGANN', 'TTGNCN', 'TTGNNA', 'TTNACN', 'TTNANA', 'TTNNCA']

# Length to extend the search window (upstream and downstream of TSS)
ten_upstream_margin = 50
ten_downstream_margin = 10

# Where to search for -35 promoters relative to -10 promoters
min_len_gap = 15
max_len_gap = 20

def search_for_promoters(promoter_list, promoter_type, TSS, window_start, window_end, ref_seq, strand):
    unique_search_results = set()
    if strand == 'F':
        seq_to_search = ref_seq[window_start:window_end]
    elif strand == 'R':
        seq_to_search = ref_seq[window_start:window_end].reverse_complement()
    else: print('strand error')
    for promoter in promoter_list:
        search_results = nt_search(str(seq_to_search), promoter)
        num_promoters_found = len(search_results) - 1
        if num_promoters_found > 0:
            for pos in search_results[1:]:
                if strand == 'F':
                    promoter_pos = window_start + pos
                    distance_to_TSS = promoter_pos + 1 - TSS
                    if promoter_type == -10:
                        promoter_seq = str(ref_seq[(promoter_pos - 3):(promoter_pos + 6)])
                    elif promoter_type == -35:
                        promoter_seq = str(ref_seq[promoter_pos:(promoter_pos + 6)])
                elif strand == 'R':
                    promoter_pos = window_end - pos
                    distance_to_TSS = TSS - promoter_pos
                    if promoter_type == -10:
                        promoter_seq = str(ref_seq[(promoter_pos - 7):(promoter_pos + 3)].reverse_complement())
                    elif promoter_type == -35:
                        promoter_seq = str(ref_seq[(promoter_pos - 6):promoter_pos].reverse_complement())
                else: print('strand error')
                unique_search_results.add((promoter_pos, distance_to_TSS, promoter_seq))
    return unique_search_results

def record_thirtyfive_ten_promoters(ref_seq, TSS, strand, upstream_margin, downstream_margin, min_len_gap, max_len_gap):
    full_results = []
    if strand == 'F':
        # Subtract one from TSS to convert to pythonic indexing.
            ### Not necessary for end position because pythonic range is 'up to but not including'
            # Define the window to look for -10 promoters
        ten_window_start_pos = TSS - 1 - ten_upstream_margin
        ten_window_end_pos = TSS + ten_downstream_margin
        max_search_win = ref_seq[(TSS - 1 - (upstream_margin + max_len_gap + 6)):(TSS + downstream_margin)]
    elif strand == 'R':
        ten_window_start_pos = TSS - 1 - ten_downstream_margin
        ten_window_end_pos = TSS + ten_upstream_margin
        max_search_win = ref_seq[(TSS - 1 - downstream_margin):(TSS + (upstream_margin + max_len_gap + 6))].reverse_complement()
    else: print('strand error')
    ten_promoter_results = search_for_promoters(ten_promoters, -10, TSS, ten_window_start_pos, ten_window_end_pos, ref_sequence, strand)
    print(ten_promoter_results)
    exists_thirtyfive_promoter_results = False
    if len(ten_promoter_results) > 0:
        for promoter_record in ten_promoter_results:
            ten_promoter_pos = promoter_record[0]
            ten_promoter_distance_to_TSS = promoter_record[1]
            ten_promoter_seq = promoter_record[2]
            if strand == 'F':
                thirtyfive_win_start = ten_promoter_pos - (max_len_gap + 6)
                thirtyfive_win_end = ten_promoter_pos - min_len_gap
            elif strand == 'R':
                thirtyfive_win_start = ten_promoter_pos + min_len_gap
                thirtyfive_win_end = the_promoter_pos + (max_len_gap + 6)
            else: print('strand error')
            thirtyfive_promoter_results = search_for_promoters(thirtyfive_promoters, -35, TSS, thirtyfive_win_start, thirtyfive_win_end, ref_sequence, strand)
            if len(thirtyfive_promoter_results) > 0:
                # -35 promoters boolean switched to true if any of the -10 promoters have -35 promoters
                exists_thirtyfive_promoter_results = True
                if ten_promoter_seq.startswith('TG'):
                    full_results.append((str(max_search_win), ten_promoter_pos, ten_promoter_distance_to_TSS, ten_promoter_seq, 'extended', thirtyfive_promoter_results))
                else:
                    full_results.append((str(max_search_win), ten_promoter_pos, ten_promoter_distance_to_TSS, ten_promoter_seq, ' ', thirtyfive_promoter_results))
            else:
                if ten_promoter_seq.startswith('TG'):
                    full_results.append((str(max_search_win), ten_promoter_pos, ten_promoter_distance_to_TSS, ten_promoter_seq, 'extended'))
                else:
                    full_results.append((str(max_search_win), ten_promoter_pos, ten_promoter_distance_to_TSS, ten_promoter_seq, ' '))
    return (full_results, exists_thirtyfive_promoter_results)

with open('test.csv', 'w') as outfile:
    writer = csv.writer(outfile, delimiter = '\t')
    writer.writerow(['uid', 'threshold', 'genome', 'feature name', 'strand', 'TSS', 'search results'])
    for entry in TSS_list:
        uid = entry[0]
        threshold = entry[1]
        genome = entry[2]
        feature_name = entry[3]
        strand = entry[4]
        TSS = entry[5]
        full_results = record_thirtyfive_ten_promoters(ref_sequence, TSS, strand, ten_upstream_margin, ten_downstream_margin, min_len_gap, max_len_gap)
        promoter_data = full_results[0]
        thirtyfive_promoter_presense = full_results[1]
        if len(promoter_data) > 0:
                # Print to results if there are any -35 promoters (which were only searched for if there were -10 promoters)
                ### or if -10 promoter is marked 'extended'
                if thirtyfive_promoter_presense is True or any('extended' in result[3] for result in promoter_data):
                    writer.writerow([uid, threshold, genome, feature_name, strand, TSS, promoter_data])