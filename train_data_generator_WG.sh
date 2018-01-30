#!/bin/bash
if [ $# -ne 5 ];
then
echo "usage: "$(basename $0) "[bam-file]" "[ref-file]" "[vcf-file]" "[output_dir]" "[max_threads]"
exit
fi

bam_file=$1
ref_file=$2
vcf_file=$3
output_dir=$4
max_threads=$5
now=$(date '+%d%m%Y%H%M%S')
output_dir=$output_dir'run-'$now
mkdir $output_dir
echo $output_dir
mkdir tmp
chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
# chrs=(3)
for i in ${chrs[@]};
    do
    current_output_dir=$output_dir/'chr'$i/
    mkdir $current_output_dir
    echo "Starting chr" 'chr'$i
    python3 train_data_generator.py --bam $bam_file --ref $ref_file --vcf $vcf_file --chromosome_name 'chr'$i --output_dir $current_output_dir --max_threads $max_threads --parallel True 2>tmp/progress-$i.txt &
    done
