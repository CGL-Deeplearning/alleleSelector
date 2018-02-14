from modules.IntervalTree import IntervalTree
from modules.TsvHandler import TsvHandler
from time import time

regions_bed_path = "/Users/saureous/data/ConfidentRegions.hg19.chr1.bed"
alleles_bed_path = ""

tsv_handler = TsvHandler(regions_bed_path)

def test_min_bucket_parameter():
    """
    iterate through a range of bucket sizes, build a tree, access every element in the tree, printing the time cost for
    building and accessing at each bucket size
    :return:
    """
    time_start = time()

    intervals = tsv_handler.get_bed_intervals()

    time0 = time()
    print("intervals loaded: ", time0-time_start)

    bucket_sizes = [1,2,4,8,16,32,64,128]

    for size in bucket_sizes:
        time1 = time()
        print(size)

        interval_tree = IntervalTree(intervals, min_bucket=size)

        time2 = time()
        print("tree built:       ", time2 - time1)

        for interval in intervals:
            start = interval[0]

            matches = interval_tree.find(start, start)

            if len(matches) == 0:
                print("WARNING no matches:", start)

        time3 = time()
        print("tree accessed:    ", time3 - time2)


def test_interval_subsetting():
    """
    use TSV reader as it would be used for a BAM, build a tree, access every element in the tree, print times
    :return:
    """
    endpoints = list(range(0, 248000000, 200000))

    for i in range(len(endpoints)-1):
        start = endpoints[i]
        stop = endpoints[i+1]

        time0 = time()

        # collect intervals from BED in illumina PG standards and convert to intervals that make sense: 0-based, closed
        bed_intervals = tsv_handler.get_subset_of_bed_intervals(start=start,
                                                                stop=stop,
                                                                universal_offset=-1,
                                                                start_offset=1)

        time1 = time()
        print("intervals loaded: ", time1-time0)

        interval_tree = IntervalTree(bed_intervals)

        time2 = time()
        print("tree built:       ", time2-time1)

        for interval in bed_intervals:
            bed_start = interval[0]
            length = interval[1]-interval[0]

            matches = interval_tree.find(bed_start, bed_start)

            if len(matches) == 0:
                print("WARNING no matches:", interval)
            if length == 0:
                print("WARNING 0 length interval:", interval)
            if length < 0:
                print("WARNING negative length interval:", interval)

        time3 = time()
        print("tree accessed:    ", time3-time2)

test_interval_subsetting()