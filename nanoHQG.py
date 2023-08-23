import argparse
from utilis import *
import os
import re


def main():
    parser = argparse.ArgumentParser(description="nanoHQG Pipeline")
    parser.add_argument("-l", "--long_read", required=True, help="Use long read sequencing")
    parser.add_argument("-1", "--forward_read", required=True, help="Path to the first input FastQ file")
    parser.add_argument("-2", "--reverse_read", required=True, help="Path to the second input FastQ file")
    parser.add_argument("-ref", "--reference_genome", required=True, help="Path to the reference genome file")
    parser.add_argument("-m", "--module", required=True, help="medaka module")
    parser.add_argument("-t", "--num_threads", type=int, default=8, help="Number of threads for processing")
    parser.add_argument("-o", "--output_directory", required=True, help="Output directory for results")

    args = parser.parse_args()

    long_read = args.long_read
    forward_read = args.forward_read
    reverse_read = args.reverse_read
    reference_genome = args.reference_genome
    num_threads = args.num_threads
    module = args.module
    output_dir = args.output_directory

    # Extracting sample names from input file names or paths
    sample_name = os.path.splitext(os.path.basename(long_read))[0].split('.')[0]
    lng_reads = os.path.abspath(long_read)
    fwd_reads = os.path.abspath(forward_read)
    rev_reads = os.path.abspath(reverse_read)
    os.mkdir(output_dir)
    os.chdir(output_dir)

    # Extract reference genome unmapping reads, these reads are microbes' reads
    aligner = MinimapAligner(reference_genome, lng_reads, sample_name + '.sam')
    aligner.run_alignment(num_threads)

    extract = SamtoolsProcessor(sample_name + '.sam')
    extract.extract_mapping_reads(num_threads, sample_name + '.microbe.fastq')

    microbe_fastq = os.path.abspath(sample_name + '.microbe.fastq')

    # Flye assembly in meta mode
    assembly = FlyeAssembler(microbe_fastq, sample_name + '.flye')
    assembly.run_assembly(num_threads)

    assembly_dir = sample_name + '.flye'
    os.chdir(assembly_dir)

    fly_assembly = 'assembly.fasta'

    # Medaka polish and add multi-processors module
    medaka_polish = MedakaConsensusRunner(microbe_fastq, fly_assembly)
    medaka_polish.run_mini_align('assembly', num_threads)

    contig_names = []

    with open('assembly.fasta', 'r+') as af, open('contig_names.txt', 'w+') as cn:
        while True:
            line = af.readline()
            if line.startswith(">"):
                print(line)
                contig_names.append(line.strip()[1:])
                cn.write(line.strip()[1:] + '\n')
            elif len(line) == 0:
                break

    medaka_polish.run_consensus_parallel('assembly.bam', contig_names, module)
    medaka_polish.run_stitch(fly_assembly, 'medaka.polished.fasta', num_threads)
    os.system(f"rm -rf *hdf")

    medaka_polish_assembly = 'medaka.polished.fasta'

    # polca polish

    aligner = MinimapAligner(medaka_polish_assembly, microbe_fastq, 'assembly.sam')
    aligner.run_alignment(num_threads)

    racon_polish = RaconCorrector(microbe_fastq, 'assembly.sam', medaka_polish_assembly, 'racon_polish.fasta')
    racon_polish.run_correction(num_threads)

    racon_final_path = "racon_polish.fasta"
    racon_final_mod_path = "racon_polish.mod.fasta"

    with open(racon_final_path, 'r') as input_file:
        content = input_file.read()

    modified_content = re.sub(r'[:,-]', '_', content)

    with open(racon_final_mod_path, 'w') as output_file:
        output_file.write(modified_content)

    # polca polish
    polca_polish = PolcaPolish('racon_polish.mod.fasta', fwd_reads, rev_reads)
    polca_polish.run_polca(num_threads)

    os.rename('racon_polish.mod.fasta.PolcaCorrected.fa', 'polca.polished.fasta')

    bwa_align = BwaAligner('polca.polished.fasta', fwd_reads, rev_reads)
    bwa_align.run_indexing()
    bwa_align.run_fwd_alignment(num_threads, 'fwd_align.sam')
    bwa_align.run_rev_alignment(num_threads, 'rev_align.sam')

    # Polypolish
    polypolish = Polypolish('polca.fasta', 'fwd_align.sam', 'rev_align.sam', 'fwd_filter.sam', 'rev_filter.sam',
                            'polypolish.fasta')
    polypolish.run_polypolish()

    minimap_align = MinimapAligner('polypolish.fasta', microbe_fastq, 'polypolish.sam')
    minimap_align.run_alignment(num_threads)

    samtools_process = SamtoolsProcessor('polypolish.sam')
    samtools_process.convert_sam_to_bam('polypolish.fasta', 'polypolish.bam', num_threads)
    samtools_process.sort_bam('polypolish.bam', 'polypolish.sort', num_threads)

    # semibin2 binning
    binning = SemiBin2Runner('polypolish.fasta', 'polypolish.sort', 'semibin_out')
    binning.run_semibin2()


if __name__ == "__main__":
    main()
