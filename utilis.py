import os
import multiprocessing

class MinimapAligner:
    def __init__(self, ref, fastq, out):
        self.ref = ref
        self.fastq = fastq
        self.out = out

    def run_alignment(self, num_threads):
        command = f"minimap2 -ax map-ont -t {num_threads} {self.ref} {self.fastq} > {self.out}"
        os.system(command)


class SamtoolsProcessor:
    def __init__(self, input_sam):
        self.input_sam = input_sam

    def convert_sam_to_bam(self, input_fasta, output_bam, num_threads):
        command = f"samtools view -@ {num_threads} -bt {input_fasta} -o {output_bam} {self.input_sam}"
        os.system(command)

    @staticmethod
    def sort_bam(input_bam, output_sorted_bam, num_threads):
        command = f"samtools sort -@ {num_threads} {input_bam} -o {output_sorted_bam}"
        os.system(command)

    def extract_mapping_reads(self, num_threads, output_fastq):
        command = f"samtools fastq -n -f 4 -@ {num_threads} {self.input_sam} > {output_fastq}"
        os.system(command)


class FlyeAssembler:
    def __init__(self, input_fastq, output_dir):
        self.input_fastq = input_fastq
        self.output_dir = output_dir

    def run_assembly(self, num_threads):
        command = f"flye --nano-hq {self.input_fastq} -o {self.output_dir} -t {num_threads} --meta"
        os.system(command)


class RaconCorrector:
    def __init__(self, input_fastq, input_sam, input_assembly, output_filename):
        self.input_fastq = input_fastq
        self.input_sam = input_sam
        self.input_assembly = input_assembly
        self.output_filename = output_filename

    def run_correction(self, num_threads):
        command = f"racon -m 8 -x -6 -g -8 -t {num_threads} {self.input_fastq} {self.input_sam} \
            {self.input_assembly} > {self.output_filename}"
        os.system(command)


class MedakaConsensusRunner:
    def __init__(self, input_fastq, input_reference):
        self.input_fastq = input_fastq
        self.input_reference = input_reference

    def run_mini_align(self, prefix, num_threads):
        command = f"mini_align -i {self.input_fastq} -r {self.input_reference} -P -m -p {prefix} -t {num_threads}"
        os.system(command)

    @staticmethod
    def run_consensus(input_bam, output_hdf, module, region):
        command = f"medaka consensus {input_bam} {output_hdf} --model {module} --batch 200 --threads 2 \
         --region {region}"
        os.system(command)

    def run_consensus_parallel(self, input_bam, contig_names, module):
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        for contig_name in contig_names:
            pool.apply_async(self.run_consensus, args=(input_bam, contig_name + '.hdf', module, contig_name))
        pool.close()
        pool.join()

    @staticmethod
    def run_stitch(draft_assembly, polised_assembly, num_threads):
        command = f"medaka stitch *.hdf {draft_assembly} {polised_assembly} --threads {num_threads}"
        os.system(command)

    def run_medaka_consensus(self, num_threads, output_directory):
        command = f"medaka_consensus -i {self.input_fastq} -d {self.input_reference} -o {output_directory} \
        -t {num_threads} -m r941_min_high_g360"
        os.system(command)


class PolcaPolish:
    def __init__(self, input_assembly, forward_reads, reverse_reads):
        self.input_assembly = input_assembly
        self.forward_reads = forward_reads
        self.reverse_reads = reverse_reads

    def run_polca(self, num_threads):
        command = f"/home/server/MaSuRCA-4.1.0/bin/polca.sh -a {self.input_assembly} \
         -r '{self.forward_reads} {self.reverse_reads}' -t {num_threads} -m 3G"
        os.system(command)


class BwaAligner:
    def __init__(self, input_reference, input_fwd_fastq, input_rev_fastq):
        self.input_reference = input_reference
        self.input_fwd_fastq = input_fwd_fastq
        self.input_rev_fastq = input_rev_fastq

    def run_indexing(self):
        command = f"bwa index {self.input_reference}"
        os.system(command)

    def run_fwd_alignment(self, num_threads, fwd_sam):
        command = f"bwa mem -t {num_threads} -a {self.input_reference} {self.input_fwd_fastq} > {fwd_sam}"
        os.system(command)

    def run_rev_alignment(self, num_threads, rev_sam):
        command = f"bwa mem -t {num_threads} -a {self.input_reference} {self.input_rev_fastq} > {rev_sam}"
        os.system(command)

    def run_alignment(self, num_threads, out_sam):
        command = f"bwa mem -t {num_threads} -a {self.input_reference} {self.input_fwd_fastq} {self.input_rev_fastq} \
        > {out_sam}"
        os.system(command)


class Polypolish:
    def __init__(self, input_reference, input_alignments_1, input_alignments_2, output_filtered_1, output_filtered_2,
                 output_polished):
        self.input_reference = input_reference
        self.input_alignments_1 = input_alignments_1
        self.input_alignments_2 = input_alignments_2
        self.output_filtered_1 = output_filtered_1
        self.output_filtered_2 = output_filtered_2
        self.output_polished = output_polished

    def run_polypolish(self):
        # Run polypolish_insert_filter.py
        command = f"polypolish_insert_filter.py --in1 {self.input_alignments_1} --in2 {self.input_alignments_2} \
        --out1 {self.output_filtered_1} --out2 {self.output_filtered_2}"
        os.system(command)

        # Run polypolish
        command = f"polypolish {self.input_reference} {self.output_filtered_1} {self.output_filtered_2} \
         > {self.output_polished}"
        os.system(command)


class SemiBin2Runner:
    def __init__(self, input_fasta, input_bam, output_directory):
        self.input_fasta = input_fasta
        self.input_bam = input_bam
        self.output_directory = output_directory

    def run_semibin2(self):
        command = f"SemiBin2 single_easy_bin --environment human_oral -i {self.input_fasta} \
                  -b {self.input_bam} -o {self.output_directory} --sequencing-type long_read"
        os.system(command)
