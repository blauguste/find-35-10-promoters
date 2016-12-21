from Bio import SeqIO, SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import nt_search
import csv

gb_path = 'NC_004350_2.gb'
ref_record = SeqIO.read(open(gb_path), 'genbank')
ref_sequence = ref_record.seq
#ref_sequence = Seq('CCCCCCCCCCTGGGGGGGGGGGGGCATTATACCCCCCCCCCTTTTTTTTTTTGGTCAATTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAA', IUPAC.unambiguous_dna)
print(len(ref_sequence))

# List of transcription start sites. 
### Csv with [0] uid, [1] threshold, [2] genome, [3] feature name, [4] strand, [5] TSS
TSS_list = []
with open('random_input.csv', 'r') as TSS_input:
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

def search_for_substring(promoter_list, seq_to_search):
    search_array = []
    for promoter in promoter_list:
        search_results = nt_search(str(seq_to_search), promoter)
        search_array.append(search_results)
    return search_array

with open('random_search_results.csv', 'w') as outfile:
    writer = csv.writer(outfile, delimiter = '\t')
    writer.writerow(['uid', 'threshold', 'genome', 'feature name', 'strand', 'TSS', 'sequence searched', 'promoter pos relative to TSS, -10 promoter seq', 'promoter pos relative to TSS, -35 promoter seq'])
    for entry in TSS_list:
        uid = entry[0]
        threshold = entry[1]
        genome = entry[2]
        feature_name = entry[3]
        strand = entry[4]
        TSS = entry[5]
        if strand == 'F':
            # Subtract one from TSS to convert to pythonic indexing.
            ### Not necessary for end position because pythonic range is 'up to but not including'
            # Define the window to look for -10 promoters
            ten_window_start_pos = TSS - 1 - ten_upstream_margin
            ten_window_end_pos = TSS + ten_downstream_margin
            ten_window_to_search = ref_sequence[ten_window_start_pos:ten_window_end_pos]
            print('feature name: ', feature_name)
            print('length of search window: ', len(ten_window_to_search))
            # Look for -10 promoter
            ten_promoter_results = search_for_substring(ten_promoters, ten_window_to_search)
            print(ten_promoter_results)
            ten_results_array = []
            for ten_pos_list in ten_promoter_results:
                num_ten_results = len(ten_pos_list) - 1
                if num_ten_results >= 1:
                    for i, ten_entry in enumerate(ten_pos_list, start=1):
                        if i < len(ten_pos_list):
                            # the promoter position in pythonic index
                            ten_promoter_pos = ten_window_start_pos + ten_pos_list[i]
                            ten_distance_to_TSS = ten_promoter_pos + 1 - TSS
                            ten_promoter_seq = ref_sequence[ten_promoter_pos:ten_promoter_pos + 7]
                            # Mark if -10 'extended' promoter
                            if ten_promoter_seq[-1] == 'G':
                                ten_results_array.append((ten_distance_to_TSS, str(ten_promoter_seq), 'extended'))
                            else:
                                ten_results_array.append((ten_distance_to_TSS, str(ten_promoter_seq), ' '))
                            # Define the window to look for -35 promoters: -15 to -20 bp from -10 promoter 
                            thirtyfive_win_start = ten_promoter_pos - 26
                            thirtyfive_win_end = ten_promoter_pos - 15
                            thirtyfive_win_to_search = ref_sequence[thirtyfive_win_start:thirtyfive_win_end]
                            print('length of search window: ', len(thirtyfive_win_to_search))
                            print(thirtyfive_win_to_search)
                            # Look for the -35 promoter 
                            thirtyfive_promoter_results = search_for_substring(thirtyfive_promoters, thirtyfive_win_to_search)
                            thirtyfive_results_array = []
                            for thirtyfive_pos_list in thirtyfive_promoter_results:
                                num_thirtyfive_results = len(thirtyfive_pos_list) - 1
                                if num_thirtyfive_results >= 1:
                                    for j, thirtyfive_entry in enumerate(thirtyfive_pos_list, start=1):
                                        if j < len(thirtyfive_pos_list):
                                            thirtyfive_promoter_pos = thirtyfive_win_start + thirtyfive_pos_list[j]
                                            thirtyfive_distance_to_TSS = thirtyfive_promoter_pos + 1 - TSS
                                            thirtyfive_promoter_seq = ref_sequence[thirtyfive_promoter_pos:thirtyfive_promoter_pos + 6]
                                            thirtyfive_results_array.append((thirtyfive_distance_to_TSS, str(thirtyfive_promoter_seq)))
            max_search_win = ref_sequence[(TSS - 1 - 76):(TSS + 10)]
            print('full_search_window: ', str(max_search_win))
            print('length_total_search_window: ', len(max_search_win))
            # Print to results if there are any -35 promoters (which were only searched for if there were -10 promoters)
            if len(ten_results_array) > 0:
                unique_ten_results = set(ten_results_array)
                if len(thirtyfive_results_array) > 0:
                    unique_thirtyfive_results = set(thirtyfive_results_array)
                    writer.writerow([uid, threshold, genome, feature_name, strand, TSS, str(max_search_win), unique_ten_results, unique_thirtyfive_results])
                elif any('extended' in results[2] for results in unique_ten_results):
                    writer.writerow([uid, threshold, genome, feature_name, strand, TSS, str(max_search_win), unique_ten_results])
        elif strand == 'R':
            # Subtract one from TSS to convert to pythonic indexing.
            ### Not necessary for end position because pythonic range is 'up to but not including'
            ten_window_start_pos = TSS - 1 - ten_downstream_margin
            ten_window_end_pos = TSS + ten_upstream_margin
            ten_window_to_search = ref_sequence[ten_window_start_pos:ten_window_end_pos].reverse_complement()
            print(ten_window_to_search)
            print('feature name: ', feature_name)
            print('length of search window: ', len(ten_window_to_search))
            ten_promoter_results = search_for_substring(ten_promoters, ten_window_to_search)
            print(ten_promoter_results)
            ten_results_array = []
            for ten_pos_list in ten_promoter_results:
                num_ten_results = len(ten_pos_list) - 1
                if num_ten_results >= 1:
                    for i, ten_entry in enumerate(ten_pos_list, start=1):
                        if i < len(ten_pos_list):
                            # the promoter position via pythonic counting
                            ten_promoter_pos = ten_window_end_pos - ten_pos_list[i]
                            ten_distance_to_TSS = TSS - ten_promoter_pos
                            ten_promoter_seq = ref_sequence[(ten_promoter_pos - 7):ten_promoter_pos].reverse_complement()
                            # Mark if -10 'extended' promoter
                            if ten_promoter_seq[-1] == 'G':
                                ten_results_array.append((ten_distance_to_TSS, str(ten_promoter_seq), 'extended'))
                            else:
                                ten_results_array.append((ten_distance_to_TSS, str(ten_promoter_seq), ' '))
                            # Define the window to look for -35 promoters: -15 to -20 bp from -10 promoter 
                            thirtyfive_win_start = ten_promoter_pos + 15
                            thirtyfive_win_end = ten_promoter_pos + 26
                            thirtyfive_win_to_search = ref_sequence[thirtyfive_win_start:thirtyfive_win_end].reverse_complement()
                            print('length of search window: ', len(thirtyfive_win_to_search))
                            print('thirtyfive search window: ', thirtyfive_win_to_search.reverse_complement())
                            # Look for the -35 promoter 
                            thirtyfive_promoter_results = search_for_substring(thirtyfive_promoters, thirtyfive_win_to_search)
                            print(thirtyfive_promoter_results)
                            thirtyfive_results_array = []
                            for thirtyfive_pos_list in thirtyfive_promoter_results:
                                num_thirtyfive_results = len(thirtyfive_pos_list) - 1
                                if num_thirtyfive_results >= 1:
                                    for j, thirtyfive_entry in enumerate(thirtyfive_pos_list, start=1):
                                        if j < len(thirtyfive_pos_list):
                                            thirtyfive_promoter_pos = thirtyfive_win_end - thirtyfive_pos_list[j]
                                            thirtyfive_distance_to_TSS = TSS - thirtyfive_promoter_pos
                                            thirtyfive_promoter_seq = ref_sequence[(thirtyfive_promoter_pos - 6):thirtyfive_promoter_pos].reverse_complement()
                                            thirtyfive_results_array.append((thirtyfive_distance_to_TSS, str(thirtyfive_promoter_seq)))
            max_search_win = ref_sequence[(TSS - 1 - 10):(TSS + 76)]
            print('max_search_win: ', max_search_win.reverse_complement())
            print('length_total_search_window: ', len(max_search_win))
            if len(ten_results_array) > 0:
                unique_ten_results = set(ten_results_array)
                if len(thirtyfive_results_array) > 0:
                    unique_thirtyfive_results = set(thirtyfive_results_array)
                    writer.writerow([uid, threshold, genome, feature_name, strand, TSS, str(max_search_win.reverse_complement()), unique_ten_results, unique_thirtyfive_results])
                elif any('extended' in results[2] for results in unique_ten_results):
                    writer.writerow([uid, threshold, genome, feature_name, strand, TSS, str(max_search_win.reverse_complement()), unique_ten_results])

