version 1.0

import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV
import "../../../tasks/Utility/Finalize.wdl" as FF

task RunMinimap2 {
    input {
        File ccs_fastq
        File reference_fasta
        String sample_name
        String reference_name
        String preset = "map-hifi"
        Int threads = 16
        String memory = "32G"
        String disk_size = "100G"
    }

    command <<<
        set -euo pipefail

        minimap2 --version
        samtools --version

        minimap2 -ax ~{preset} \
            -t ~{threads} \
            ~{reference_fasta} \
            ~{ccs_fastq} | \
        samtools sort -@ ~{threads} -O BAM -o ~{sample_name}.~{reference_name}.sorted.bam
        
        samtools index -@ ~{threads} ~{sample_name}.~{reference_name}.sorted.bam
    >>>

    output {
        File aligned_bam = "~{sample_name}.~{reference_name}.sorted.bam"
        File aligned_bai = "~{sample_name}.~{reference_name}.sorted.bam.bai"
    }

    runtime {
        docker: "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0"
        memory: memory
        cpu: threads
        disks: "local-disk " + disk_size + " SSD"
    }
}

workflow Minimap2AlignmentWithMetrics {
    input {
        File ccs_fastq
        File reference_fasta
        String reference_name
        String sample_name
        String preset = "map-hifi"
        String outdir
        File? bed_to_compute_coverage
    }

    call RunMinimap2 {
        input:
            ccs_fastq = ccs_fastq,
            reference_fasta = reference_fasta,
            sample_name = sample_name,
            reference_name = reference_name,
            preset = preset
    }

    call COV.SampleLevelAlignedMetrics as coverage {
        input:
            aligned_bam = RunMinimap2.aligned_bam,
            aligned_bai = RunMinimap2.aligned_bai,
            ref_fasta = reference_fasta,
            bed_to_compute_coverage = bed_to_compute_coverage
    }

    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeBam { 
        input: 
            outdir = dir, 
            file = RunMinimap2.aligned_bam, 
            name = "~{sample_name}.~{reference_name}.bam" 
    }
    call FF.FinalizeToFile as FinalizeBai { 
        input: 
            outdir = dir, 
            file = RunMinimap2.aligned_bai, 
            name = "~{sample_name}.~{reference_name}.bam.bai" 
    }

    output {
        File ~{reference_name}_aligned_bam = FinalizeBam.gcs_path
        File ~{reference_name}_aligned_bai = FinalizeBai.gcs_path

        Float ~{reference_name}_aligned_num_reads = coverage.aligned_num_reads
        Float ~{reference_name}_aligned_num_bases = coverage.aligned_num_bases
        Float ~{reference_name}_aligned_frac_bases = coverage.aligned_frac_bases
        Float ~{reference_name}_aligned_est_fold_cov = coverage.aligned_est_fold_cov
        Float ~{reference_name}_aligned_read_length_mean = coverage.aligned_read_length_mean
        Float ~{reference_name}_aligned_read_length_median = coverage.aligned_read_length_median
        Float ~{reference_name}_aligned_read_length_stdev = coverage.aligned_read_length_stdev
        Float ~{reference_name}_aligned_read_length_n50 = coverage.aligned_read_length_N50
        Float ~{reference_name}_average_identity = coverage.average_identity
        Float ~{reference_name}_median_identity = coverage.median_identity
    }
}

