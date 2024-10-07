version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Assembly/Hifiasm.wdl" as HA
import "../../../tasks/QC/Quast.wdl" as QuastEval
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBAssembleWithHifiasm {

    meta {
        description: "A workflow that performs single sample genome assembly on PacBio HiFi reads from one or more SMRT cells. The multiple SMRT cells data are merged prior to assembly."
    }
    parameter_meta {
        ccs_fqs:            "GCS path to CCS fastq files"

        participant_name:   "name of the participant from whom these samples were obtained"
        prefix:             "prefix for output files"

        ref_fasta_for_eval: "Reference Fasta used for evaluating "
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        Array[File] ccs_fqs

        String participant_name
        String prefix

        File? ref_fasta_for_eval

        String gcs_out_root_dir
    }

    #########################################################################################
    if (length(ccs_fqs) > 1) {
        call Utils.MergeFastqs as MergeAllFastqs { input: fastqs = ccs_fqs }
    }
    File ccs_fq  = select_first([ MergeAllFastqs.merged_fastq, ccs_fqs[0] ])

    call HA.Hifiasm {
        input:
            reads = ccs_fq,
            prefix = prefix
    }

    #########################################################################################
    # Finalize data
    String workflow_name = "PBAssembleWithHifiasm"

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


    # assembly results themselves
    call FF.CompressAndFinalize as FinalizeHifiasmPrimaryGFA   { input: outdir = dir, file = Hifiasm.primary_gfa }
    call FF.CompressAndFinalize as FinalizeHifiasmPrimaryFA    { input: outdir = dir, file = Hifiasm.primary_tigs }

    call FF.CompressAndFinalize as FinalizeHifiasmAlternateGFA   { input: outdir = dir, file = Hifiasm.alternate_gfa }
    call FF.CompressAndFinalize as FinalizeHifiasmAlternateFA    { input: outdir = dir, file = Hifiasm.alternate_tigs }

    call FF.FinalizeAndCompress as FinalizeHifiasmHapGFAs  { input: outdir = dir, files = Hifiasm.phased_gfas, prefix = prefix + ".haploGFAs" }
    call FF.FinalizeAndCompress as FinalizeHifiasmHapFAs   { input: outdir = dir, files = Hifiasm.phased_tigs, prefix = prefix + ".haploTigs" }

    output {
        File merged_fq = finalized_merged_fq_path

        File hifiasm_primary_gfa  = FinalizeHifiasmPrimaryGFA.gcs_path
        File hifiasm_primary_tigs = FinalizeHifiasmPrimaryFA.gcs_path

        File hifiasm_haploGFAs = FinalizeHifiasmHapGFAs.gcs_path
        File hifiasm_haplotigs = FinalizeHifiasmHapFAs.gcs_path

        File hifiasm_alternate_gfa  = FinalizeHifiasmAlternateGFA.gcs_path
        File hifiasm_alternate_tigs = FinalizeHifiasmAlternateFA.gcs_path
    }
}
