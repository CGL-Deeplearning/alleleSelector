# alleleSelector

From a genome alignment, this program finds possible candidate windows for variant calling.

## Dependencies
python3.5+, pysam

#### Example run:
python3 main.py --bam [path_to_bam] --ref [path to reference fasta file] --chromosome_name chr3 --max_threads [max_number_of_threads] --test [True/False] --jsoon [True/False]

