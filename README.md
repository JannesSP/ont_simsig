# squigseg
First steps of creating a better resquiggling algorithm than the existing ones


## Basecall testdata

### Home

    "C:\Program Files\OxfordNanopore\ont-guppy\bin\guppy_basecaller.exe" -i data\simulation\ -s data\basecalls\ --disable_pings --device cuda:0 --disable_qscore_filtering --calib_detect -c "C:\Program Files\OxfordNanopore\ont-guppy\data\rna_r9.4.1_70bps_hac.cfg"

### Work

    ~/bin/ont-guppy/bin/guppy_basecaller -i data/simulation/ -s data/basecalls/ --disable_pings --device cuda:0 --disable_qscore_filtering --calib_detect -c ~/bin/ont-guppy/data/rna_r9.4.1_70bps_hac.cfg

## map basecalls

    minimap2 -ax map-ont -k 14 -t 20 data/simulation/RNA_simulation_20220704_171710/RNA_simulation_20220704_171710_batch_reference.fasta data/basecalls/RNA_simulation_20220704_171710/fastq_runid_20220704_simulation_run_id_0_0.fastq | samtools view -hbF4 | samtools sort > data/basecalls/RNA_simulation_20220704_171710/fastq.bam

    samtools index data/basecalls/RNA_simulation_20220704_171710/fastq.bam

    /home/yi98suv/anaconda3/envs/pysam/bin/python /home/yi98suv/projects/squigseg/src/count_mapping.py data/basecalls/RNA_simulation_20220704_171710/fastq.bam data/simulation/RNA_simulation_20220704_171710/RNA_simulation_20220704_171710_batch_reference.fasta

## Basecall new simulated data

    for dir in data/simulation/RNA_simulation_*; do ~/bin/ont-guppy/bin/guppy_basecaller -i data/simulation/$(basename $dir) -s data/basecalls/$(basename $dir) --disable_pings --device cuda:0 --disable_qscore_filtering --calib_detect -c ~/bin/ont-guppy/data/rna_r9.4.1_70bps_hac.cfg; done

Renaming fastqs

    for file in data/basecalls/*/*fastq.temp; do mv $file ${file%%.temp}; done

## Mapping

    for dir in data/simulation/RNA_simulation_*; do minimap2 -ax map-ont -k 14 -t 24 $dir/$(basename $dir)_batch_reference.fasta data/basecalls/$(basename $dir)/fastq_runid*.fastq | samtools view -hbF4 | samtools sort > data/basecalls/$(basename $dir)/fastq.bam && samtools index data/basecalls/$(basename $dir)/fastq.bam; done

---

    for dir in data/basecalls/*; do /home/yi98suv/anaconda3/envs/pysam/bin/python /home/yi98suv/projects/squigseg/src/count_mapping.py $dir/fastq.bam data/simulation/$(basename $dir)/$(basename $dir)_batch_reference.fasta; done