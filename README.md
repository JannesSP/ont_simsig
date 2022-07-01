# squigseg
First steps of creating a better resquiggling algorithm than the existing ones


## Basecall testdata

### Home

"C:\Program Files\OxfordNanopore\ont-guppy\bin\guppy_basecaller.exe" -i data\simulation\ -s data\simulation\ --disable_pings --device cuda:0 --disable_qscore_filtering --calib_detect -c "C:\Program Files\OxfordNanopore\ont-guppy\data\rna_r9.4.1_70bps_hac.cfg"

### Work

/data/mahlzeitlocal/guppy/ont-guppy/bin/guppy_basecaller -i data/simulation/ -s data/simulation/ --disable_pings --device cuda:0 --disable_qscore_filtering --calib_detect -c /data/mahlzeitlocal/guppy/ont-guppy/data/rna_r9.4.1_70bps_hac.cfg