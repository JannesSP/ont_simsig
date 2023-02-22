#Simulate simple Oxford Nanopore Technologies RNA sequencing signals

## Basecall simulated data
for dir in data/RNA_simulation_*/simulation; do ~/bin/ont-guppy/bin/guppy_basecaller -i $dir -s $(basename $dir)/basecalls/ --disable_pings --device cuda:0 --disable_qscore_filtering --calib_detect -c ~/bin/ont-guppy/data/rna_r9.4.1_70bps_hac.cfg; done

#Renaming fastqs if needed
for file in data/*/basecalls/*fastq.temp; do mv $file ${file%%.temp}; done

## Mapping
for dir in data/RNA_simulation_*/simulation; do minimap2 -ax map-ont -k 14 -t 24 $dir/$(basename $(dirname $dir))_batch_reference.fasta $(dirname $dir)/basecalls/*.fastq | samtools view -hbF4 | samtools sort > $(dirname $dir)/basecalls/fastq.bam && samtools index $(dirname $dir)/basecalls/fastq.bam; done

#---
for dir in data/RNA*; do ~/anaconda3/envs/pysam/bin/python ~/projects/squigseg/src/count_mapping.py $dir/basecalls/fastq.bam $dir/simulation/$(basename $dir)_batch_reference.fasta; done

#The first three simulations using the set_segment_length of 50 and stdev_scale of 0, 0.2 and 0.4 did not work!
#Guppy was able to basecall, but minimap2 was not able to map these reads - probably too many errors?

## running nanopolish eventalign on simulated the data
for dir in data/RNA*; do nanopolish index -d $dir/simulation -s $dir/basecalls/sequencing_summary.txt $dir/basecalls/*fastq; done
for dir in data/RNA*; do mkdir $dir/nanopolish; done
for dir in data/RNA*; do nanopolish eventalign --reads $dir/basecalls/*fastq --bam $dir/basecalls/fastq.bam --genome $dir/simulation/*_reference.fasta --summary=$dir/nanopolish/eventalign_summary.csv --scale-events --signal-index -t 24 --progress > $dir/nanopolish/eventalign_result.csv; done

## extract nanopolish segmentation
for dir in data/RNA*; do python src/extract_raw_nanopolish_seg.py $dir/nanopolish/eventalign_summary.csv $dir/nanopolish/eventalign_result.csv $dir/nanopolish/segmentation.hdf5; done

# ont_simsig.py
# Simulate 1 read of a random reference with 2000 bases
python src/ont_simsig.py 1 -rl 2000 data/simulation/ --fullRef