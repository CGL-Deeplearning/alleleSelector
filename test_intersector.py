from modules.Intersector import *

query_intervals_file_path = "/Users/saureous/data/intersector/a.bed"
reference_intervals_file_path = "/Users/saureous/data/intersector/b.bed"

intervals = intersect(bed_file_path_query=query_intervals_file_path, bed_file_path_reference=reference_intervals_file_path)

for chr in intervals:
    print(chr, intervals[chr])

# query_intervals_file_path = "/Users/saureous/data/candidate_test/test/Label_chr1_2000000_2200000_2.bed"
# reference_intervals_file_path = "/Users/saureous/data/ConfidentRegions.hg19.chr1.bed"
#
# intersect(bed_file_path_query=query_intervals_file_path, bed_file_path_reference=reference_intervals_file_path)

