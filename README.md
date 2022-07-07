# squigseg
First steps of creating a better resquiggling algorithm than the existing ones

## Basecall simulated data

    for dir in data/RNA_simulation_*/simulation; do ~/bin/ont-guppy/bin/guppy_basecaller -i $dir -s $(basename $dir)/basecalls/ --disable_pings --device cuda:0 --disable_qscore_filtering --calib_detect -c ~/bin/ont-guppy/data/rna_r9.4.1_70bps_hac.cfg; done

Renaming fastqs

    for file in data/*/basecalls/*fastq.temp; do mv $file ${file%%.temp}; done

## Mapping

    for dir in data/RNA_simulation_*/simulation; do minimap2 -ax map-ont -k 14 -t 24 $dir/$(basename $(dirname $dir))_batch_reference.fasta $(dirname $dir)/basecalls/*.fastq | samtools view -hbF4 | samtools sort > $(dirname $dir)/basecalls/fastq.bam && samtools index $(dirname $dir)/basecalls/fastq.bam; done

---

    for dir in data/RNA*; do /home/yi98suv/anaconda3/envs/pysam/bin/python /home/yi98suv/projects/squigseg/src/count_mapping.py $dir/basecalls/fastq.bam $dir/simulation/$(basename $dir)_batch_reference.fasta; done