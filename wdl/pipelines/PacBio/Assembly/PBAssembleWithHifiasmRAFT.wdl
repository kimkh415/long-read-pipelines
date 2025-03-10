version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Assembly/Hifiasm.wdl" as HA
import "../../../tasks/Assembly/HifiasmRAFT.wdl" as RAFT
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBAssembleWithHifiasm {

    meta {
        description: "A workflow that performs single sample genome assembly on PacBio HiFi reads from one or more SMRT cells. The multiple SMRT cells data are merged prior to assembly."
    }
    parameter_meta {
        ccs_fqs:            "GCS path to CCS fastq files"
        genome_length:      "Estimated length of genome in bases"
        est_coverage:       "estimated fold coverage of the reads"
        error_correct_read_fa: "from hifiasm error correct reads"
        overlaps_paf:       "from hifiasm all to all alignment overlaps"
        prefix:             "prefix for output files"
        raft_disk_size:     "Disk size in GB [300]"
        raft_memory:        "Memory in GB [200]"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        Array[File] ccs_fqs
        String genome_length = "2922918302"
        Int? est_coverage
        File? error_correct_read_fa
        File? overlaps_paf

        String prefix
        Int raft_disk_size = 300
        Int raft_memory = 200

        String gcs_out_root_dir
    }

    #########################################################################################
    if (length(ccs_fqs) > 1) {
        call Utils.MergeFastqs as MergeAllFastqs { input: fastqs = ccs_fqs }
    }
    File ccs_fq  = select_first([ MergeAllFastqs.merged_fastq, ccs_fqs[0] ])

    call RAFT.FragmentReadsRAFT {
        input:
            reads = ccs_fq,
            genome_length = genome_length,
            est_coverage = est_coverage,
            error_correct_read_fa = error_correct_read_fa,
            overlaps_paf = overlaps_paf,
            raft_disk_size = raft_disk_size,
            raft_memory = raft_memory
    }

    call HA.Hifiasm {
        input:
            reads = FragmentReadsRAFT.fragmented_reads,
            prefix = prefix
    }

    #########################################################################################
    # Finalize data
    String workflow_name = "PBAssembleWithHifiasmRAFT"

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name + "/~{prefix}"
    String dir = outdir + "/assembly"

    # merged FASTQ
    String dummy = basename(ccs_fq)
    String dummy_b = sub(dummy, ".gz$", "")
    if (dummy != dummy_b) {
        call FF.FinalizeToFile as FinalizeMergedFQ { input: outdir = dir, file = ccs_fq, name = prefix + ".fq.gz" }
    }
    if (dummy == dummy_b) {
        call FF.CompressAndFinalize as CompressAndFinalizeMergedFQ { input: outdir = dir, file = ccs_fq, name = prefix + ".fq.gz" }
    }
    String finalized_merged_fq_path = select_first([FinalizeMergedFQ.gcs_path, CompressAndFinalizeMergedFQ.gcs_path])

    # RAFT results
    call FF.CompressAndFinalize as FinalizeRaftECReads  { input: outdir = outdir + "/RAFT", file = FragmentReadsRAFT.ec_reads }
    call FF.CompressAndFinalize as FinalizeRaftOverlaps  { input: outdir = outdir + "/RAFT", file = FragmentReadsRAFT.overlaps }

    # assembly results themselves
    call FF.CompressAndFinalize as FinalizeHifiasmPrimaryGFA   { input: outdir = dir, file = Hifiasm.primary_gfa }
    call FF.CompressAndFinalize as FinalizeHifiasmPrimaryFA    { input: outdir = dir, file = Hifiasm.primary_tigs }

    call FF.CompressAndFinalize as FinalizeHifiasmAlternateGFA   { input: outdir = dir, file = Hifiasm.alternate_gfa }
    call FF.CompressAndFinalize as FinalizeHifiasmAlternateFA    { input: outdir = dir, file = Hifiasm.alternate_tigs }

    call FF.FinalizeAndCompress as FinalizeHifiasmHapGFAs  { input: outdir = dir, files = Hifiasm.phased_gfas, prefix = prefix + ".haploGFAs" }
    call FF.FinalizeAndCompress as FinalizeHifiasmHapFAs   { input: outdir = dir, files = Hifiasm.phased_tigs, prefix = prefix + ".haploTigs" }

    call FF.CompressAndFinalize as FinalizeECReads  { input: outdir = dir, file = FragmentReadsRAFT.ec_reads}
    call FF.CompressAndFinalize as FinalizeOverlaps { input: outdir = dir, file = FragmentReadsRAFT.overlaps}

    output {
        File merged_fq = finalized_merged_fq_path

        Int hifiasm_raft_est_cov = FragmentReadsRAFT.coverage
        File hifiasm_raft_ec_reads = FinalizeRaftECReads.gcs_path
        File hifiasm_raft_overlaps = FinalizeRaftOverlaps.gcs_path

        File hifiasm_raft_primary_gfa  = FinalizeHifiasmPrimaryGFA.gcs_path
        File hifiasm_raft_primary_tigs = FinalizeHifiasmPrimaryFA.gcs_path

        File hifiasm_raft_haploGFAs = FinalizeHifiasmHapGFAs.gcs_path
        File hifiasm_raft_haplotigs = FinalizeHifiasmHapFAs.gcs_path

        File hifiasm_raft_alternate_gfa  = FinalizeHifiasmAlternateGFA.gcs_path
        File hifiasm_raft_alternate_tigs = FinalizeHifiasmAlternateFA.gcs_path
    }
}
