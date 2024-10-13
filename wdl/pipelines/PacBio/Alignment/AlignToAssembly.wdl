version 1.0

import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV
import "../../../tasks/Alignment/AlignReads.wdl" as AR
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow Minimap2AlignmentWithMetrics {
    input {
        File ccs_fq
        File reference_fasta_hap1
        File reference_fasta_hap2
        String sample_name
        String lib_name = "RevioCCS"
        String preset = "map-hifi"
        String outdir
        File? bed_to_compute_coverage
    }

    call AR.Minimap2 as AlignHap1 {
        input:
            reads         = ccs_fq,
            ref_fasta   = reference_fasta_hap1,
            prefix = sample_name,
            library     = lib_name,
            map_preset  = preset,
            RG          = "@RG\tID:" + sample_name + "\tPL:PacBio\tLB:" + lib_name + "\tSM:" + sample_name
    }

    call PB.Align as AlignHap2 {
        input:
            reads         = ccs_fq,
            ref_fasta   = reference_fasta_hap2,
            prefix = sample_name,
            library     = lib_name,
            map_preset  = preset,
            RG          = "@RG\tID:" + sample_name + "\tPL:PacBio\tLB:" + lib_name + "\tSM:" + sample_name
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

