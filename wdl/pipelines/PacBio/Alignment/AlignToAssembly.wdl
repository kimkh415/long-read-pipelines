version 1.0

import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV
import "../../../tasks/Alignment/AlignReads.wdl" as AR
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow Pbmm2AlignmentWithMetrics {
    input {
        File ccs_fq
        File reference_fasta_hap1
        File reference_fasta_hap2
        String sample_name
        String lib_name = "RevioCCS"
        String preset = "CCS"
        String outdir
        File? bed_to_compute_coverage
    }

    Array[File] reference_fastas = [reference_fasta_hap1, reference_fasta_hap2]

    Boolean any_gzipped = length(select_all([ 
        if sub(reference_fasta_hap1, ".*\\.", "") == "gz" then true else false,
        if sub(reference_fasta_hap2, ".*\\.", "") == "gz" then true else false
    ])) > 0

    if (any_gzipped) {
        call UnzipReference {
            input:
                gzipped_fastas = reference_fastas
        }
    }

    File unzipped_hap1_fasta = select_first([UnzipReference.unzipped_fastas[0], reference_fasta_hap1])
    File unzipped_hap2_fasta = select_first([UnzipReference.unzipped_fastas[1], reference_fasta_hap2])

    call PB.Align as AlignHap1 {
        input:
            bam         = ccs_fq,
            ref_fasta   = unzipped_hap1_fasta,
            sample_name = sample_name,
            library     = lib_name,
            map_preset  = preset,
            drop_per_base_N_pulse_tags = true
    }

    call PB.Align as AlignHap2 {
        input:
            bam         = ccs_fq,
            ref_fasta   = unzipped_hap2_fasta,
            sample_name = sample_name,
            library     = lib_name,
            map_preset  = preset,
            drop_per_base_N_pulse_tags = true
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

task UnzipReference {  
    input {  
        Array[File] gzipped_fastas  
    }  

    command <<<  
        set -e  
        mkdir unzipped_files  
        for gzipped_file in ~{sep=' ' gzipped_fastas}; do  
            output_filename=$(basename "$gzipped_file" .gz)  
            gunzip -c "$gzipped_file" > "unzipped_files/$output_filename"  
            echo "unzipped_files/$output_filename" >> unzipped_file_list.txt  
        done  
    >>>  

    output {  
        Array[File] unzipped_fastas = read_lines("unzipped_file_list.txt")  
    }  

    runtime {  
        cpu: 2  
        memory: "20 GiB"  
        disks: "local-disk 20 HDD"  
        docker: "ubuntu:20.04"  
    }  
}  
