# nanoHQG
Pipline for assembling high quality draft genomes based on low coverage ONT reads and short reads
## installation
```
git clone https://github.com/bikmi/nanoHQG
cd nanoHQG
conda env create -f nanoHQG.yaml
or
mamba create -f nanoHQG.yaml
```

## download database
human genome is required for decontamination
you can download from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/

## run 
python nanoHQG.py -l long_reads -1 fwd_short_reads -2 rev_short_reads -ref human_genome -m module -t 30 -o out_dir

usage: nanoHQG.py [-h] -l LONG_READ -1 FORWARD_READ -2 REVERSE_READ -ref
                  REFERENCE_GENOME -m MODULE [-t NUM_THREADS] -o
                  OUTPUT_DIRECTORY

nanoHQG Pipeline

optional arguments:
  -h, --help            show this help message and exit
  -l LONG_READ, --long_read LONG_READ
                        Use long read sequencing
  -1 FORWARD_READ, --forward_read FORWARD_READ
                        Path to the first input FastQ file
  -2 REVERSE_READ, --reverse_read REVERSE_READ
                        Path to the second input FastQ file
  -ref REFERENCE_GENOME, --reference_genome REFERENCE_GENOME
                        Path to the reference genome file
  -m MODULE, --module MODULE
                        medaka module
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads for processing
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Output directory for results
