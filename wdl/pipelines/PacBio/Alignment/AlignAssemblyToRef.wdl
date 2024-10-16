version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow AlignAssemblyToRef {
    input {
        File assembly_fasta1
        File assembly_fasta2
        File reference_fasta
        String sample_name
        String ref_name
        String preset = "asm5"
        String outdir
    }

    call RunMinimap2 as AlignHap1 {
        input:
            assembly_fasta = assembly_fasta1,
            reference_fasta = reference_fasta,
            sample_name = sample_name,
            haplotype = "hap1",
            preset = preset
    }

    call RunMinimap2 as AlignHap2 {
        input:
            assembly_fasta = assembly_fasta2,
            reference_fasta = reference_fasta,
            sample_name = sample_name,
            haplotype = "hap2",
            preset = preset
    }

    String dir = outdir + "/contig_alignments"

    call FF.FinalizeToFile as FinalizeBamHap1 {
        input:
            outdir = dir,
            file = AlignHap1.aligned_bam,
            name = "~{sample_name}.hap1.~{ref_name}.bam"
    }
    call FF.FinalizeToFile as FinalizeBaiHap1 {
        input:
            outdir = dir,
            file = AlignHap1.aligned_bai,
            name = "~{sample_name}.hap1.~{ref_name}.bam.bai"
    }

    call FF.FinalizeToFile as FinalizeBamHap2 {
        input:
            outdir = dir,
            file = AlignHap2.aligned_bam,
            name = "~{sample_name}.hap2.~{ref_name}.bam"
    }
    call FF.FinalizeToFile as FinalizeBaiHap2 {
        input:
            outdir = dir,
            file = AlignHap2.aligned_bai,
            name = "~{sample_name}.hap2.~{ref_name}.bam.bai"
    }

    output {
        File hap1_ref_aligned_bam = FinalizeBamHap1.gcs_path
        File hap1_ref_aligned_bai = FinalizeBaiHap1.gcs_path
        File hap2_ref_aligned_bam = FinalizeBamHap2.gcs_path
        File hap2_ref_aligned_bai = FinalizeBaiHap2.gcs_path
    }
}


task RunMinimap2 {
    input {
        File assembly_fasta
        File reference_fasta
        String sample_name
        String haplotype
        String preset
        Int threads = 16
        Int memory = 80
        Int disk_size = 100
    }

    command <<<
        set -euo pipefail

        minimap2 --version
        samtools --version

        minimap2 -ax ~{preset} \
            -t ~{threads} \
            ~{reference_fasta} \
            ~{assembly_fasta} \
            > ~{sample_name}.~{haplotype}.sam

        samtools sort -@ ~{threads} -O BAM -o ~{sample_name}.~{haplotype}.sorted.bam ~{sample_name}.~{haplotype}.sam
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
        disks: "local-disk " + disk_size + " HDD"
    }
}

