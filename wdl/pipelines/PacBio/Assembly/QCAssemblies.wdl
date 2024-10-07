version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/QC/Quast.wdl" as QuastEval
import "../../../tasks/Utility/Finalize.wdl" as FF


workflow QCAssemblies {

    meta {
        description: "Perform Quast QC on two haplotype resolved assemblies"
    }
    parameter_meta {
        ccs_fq:            "GCS path to CCS fastq file"
        ref_fasta_for_eval: "Reference Fasta used for evaluating "
        assemblies:         "list of assemblies (e.g. hap1 hap2)"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        File ccs_fq
        String prefix
        File ref_fasta_for_eval
        Array[File] assemblies
        String gcs_out_root_dir
    }

    #########################################################################################
    call QuastEval.Quast as hap_quast {
        input:
            ref = ref_fasta_for_eval,
            fq = ccs_fq,
            assemblies = assemblies
    }

    call QuastEval.SummarizeQuastReport as hap_quast_summary {
        input: quast_report_txt = hap_quast.report_txt
    }

    #########################################################################################
    # Finalize data
    String workflow_name = "PBAssembleWithHifiasm"

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name + "/~{prefix}"
    String dir = outdir + "/assembly"

    # collect results
    call FF.FinalizeToFile as FinalizeQuastReportHtml {
        input: outdir = dir, file = hap_quast.report_html
    }
    call FF.FinalizeAndCompress as FinalizeQuastReports {
        input: outdir = dir, files = hap_quast.report_in_various_formats, prefix = prefix + ".quast_reports"
    }
    call FF.FinalizeToFile as FinalizeQuastSummaryAll {
        input: outdir = dir, file = hap_quast_summary.quast_metrics_together
    }
    scatter (report in hap_quast_summary.quast_metrics ) {
        call FF.FinalizeToFile as FinalizeQuastIndividualSummary  { input: outdir = dir, file = report }
    }

    output {
        File? quast_report_html = FinalizeQuastReportHtml.gcs_path
        File? quast_report_in_various_formats = FinalizeQuastReports.gcs_path

        File? quast_summary_on_all = FinalizeQuastSummaryAll.gcs_path

        File? quast_summary_on_H0 = FinalizeQuastIndividualSummary.gcs_path[1]
        File? quast_summary_on_H1 = FinalizeQuastIndividualSummary.gcs_path[2]
    }
}
