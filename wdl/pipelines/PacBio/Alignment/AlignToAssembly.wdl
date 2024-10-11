version 1.0

import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV
import "../../../tasks/Utility/Finalize.wdl" as FF

task RunMinimap2 {
    input {
        File ccs_fastq
        File reference_fasta
        String sample_name
        String haplotype
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
        samtools sort -@ ~{threads} -O BAM -o ~{sample_name}.~{haplotype}.sorted.bam

        samtools index -@ ~{threads} ~{sample_name}.~{haplotype}.sorted.bam
    >>>

    output {
        File aligned_bam = "~{sample_name}.~{haplotype}.sorted.bam"
        File aligned_bai = "~{sample_name}.~{haplotype}.sorted.bam.bai"
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
        File reference_fasta_hap1
        File reference_fasta_hap2
        String sample_name
        String preset = "map-hifi"
        String outdir
        File? bed_to_compute_coverage
    }

    call RunMinimap2 as AlignHap1 {
        input:
            ccs_fastq = ccs_fastq,
            reference_fasta = reference_fasta_hap1,
            sample_name = sample_name,
            haplotype = "hap1",
            preset = preset
    }

    call RunMinimap2 as AlignHap2 {
        input:
            ccs_fastq = ccs_fastq,
            reference_fasta = reference_fasta_hap2,
            sample_name = sample_name,
            haplotype = "hap2",
            preset = preset
    }

    call COV.SampleLevelAlignedMetrics as CoverageHap1 {
        input:
            aligned_bam = AlignHap1.aligned_bam,
            aligned_bai = AlignHap1.aligned_bai,
            ref_fasta = reference_fasta_hap1,
            bed_to_compute_coverage = bed_to_compute_coverage
    }

    call COV.SampleLevelAlignedMetrics as CoverageHap2 {
        input:
            aligned_bam = AlignHap2.aligned_bam,
            aligned_bai = AlignHap2.aligned_bai,
            ref_fasta = reference_fasta_hap2,
            bed_to_compute_coverage = bed_to_compute_coverage
    }

    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeBamHap1 {
        input:
            outdir = dir,
            file = AlignHap1.aligned_bam,
            name = "~{sample_name}.hap1.bam"
    }
    call FF.FinalizeToFile as FinalizeBaiHap1 {
        input:
            outdir = dir,
            file = AlignHap1.aligned_bai,
            name = "~{sample_name}.hap1.bam.bai"
    }

    call FF.FinalizeToFile as FinalizeBamHap2 {
        input:
            outdir = dir,
            file = AlignHap2.aligned_bam,
            name = "~{sample_name}.hap2.bam"
    }
    call FF.FinalizeToFile as FinalizeBaiHap2 {
        input:
            outdir = dir,
            file = AlignHap2.aligned_bai,
            name = "~{sample_name}.hap2.bam.bai"
    }

    output {
        File hap1_aligned_bam = FinalizeBamHap1.gcs_path
        File hap1_aligned_bai = FinalizeBaiHap1.gcs_path
        File hap2_aligned_bam = FinalizeBamHap2.gcs_path
        File hap2_aligned_bai = FinalizeBaiHap2.gcs_path

        Float hap1_aligned_num_reads = CoverageHap1.aligned_num_reads
        Float hap1_aligned_num_bases = CoverageHap1.aligned_num_bases
        Float hap1_aligned_frac_bases = CoverageHap1.aligned_frac_bases
        Float hap1_aligned_est_fold_cov = CoverageHap1.aligned_est_fold_cov
        Float hap1_aligned_read_length_mean = CoverageHap1.aligned_read_length_mean
        Float hap1_aligned_read_length_median = CoverageHap1.aligned_read_length_median
        Float hap1_aligned_read_length_stdev = CoverageHap1.aligned_read_length_stdev
        Float hap1_aligned_read_length_n50 = CoverageHap1.aligned_read_length_N50
        Float hap1_average_identity = CoverageHap1.average_identity
        Float hap1_median_identity = CoverageHap1.median_identity

        Float hap2_aligned_num_reads = CoverageHap2.aligned_num_reads
        Float hap2_aligned_num_bases = CoverageHap2.aligned_num_bases
        Float hap2_aligned_frac_bases = CoverageHap2.aligned_frac_bases
        Float hap2_aligned_est_fold_cov = CoverageHap2.aligned_est_fold_cov
        Float hap2_aligned_read_length_mean = CoverageHap2.aligned_read_length_mean
        Float hap2_aligned_read_length_median = CoverageHap2.aligned_read_length_median
        Float hap2_aligned_read_length_stdev = CoverageHap2.aligned_read_length_stdev
        Float hap2_aligned_read_length_n50 = CoverageHap2.aligned_read_length_N50
        Float hap2_average_identity = CoverageHap2.average_identity
        Float hap2_median_identity = CoverageHap2.median_identity
    }
}

