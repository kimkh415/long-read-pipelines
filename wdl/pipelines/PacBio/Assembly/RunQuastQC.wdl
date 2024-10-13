version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/QC/Quast.wdl" as QuastEval
import "../../../tasks/Utility/Finalize.wdl" as FF


workflow RunQuastQC {

    meta {
        description: "Perform Quast QC on a primary and two haplotype resolved assemblies"
    }
    parameter_meta {
        prefix:						 "sample name"
        ref_fasta_for_eval: "Reference Fasta used for evaluating "
        assembly_primary:   "primary assembly fa"
        assembly_hap1:         "hap1 assembly fa"
        assembly_hap2:         "hap2 assembly fa"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        String prefix
        File ref_fasta_for_eval
        File assembly_primary
        File assembly_hap1
        File assembly_hap2
        String gcs_out_root_dir
    }

    Boolean is_gzipped = sub(ref_fasta_for_eval, ".*\\.", "") == "gz"

    if (is_gzipped) {
        call UnzipReference {
            input:
                gzipped_fasta = ref_fasta_for_eval
        }
    }

    File ref_fasta = select_first([UnzipReference.unzipped_fasta, ref_fasta_for_eval])

    #########################################################################################
    call QuastEval.Quast as hap_quast {
        input:
            ref = ref_fasta,
            is_large = true,
            assemblies = [assembly_primary, assembly_hap1, assembly_hap2]
    }

    call QuastEval.SummarizeQuastReport as hap_quast_summary {
        input: quast_report_txt = hap_quast.report_txt
    }

    #########################################################################################
    # Finalize data
    String workflow_name = "RunQuastQC"

    String dir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name + "/~{prefix}"

    # collect results
    call FF.FinalizeToFile as FinalizeQuastReportHtml {
        input: outdir = dir, file = hap_quast.report_html
    }
    call FF.FinalizeToFile as FinalizeQuastReportTxt {
        input: outdir = dir, file = hap_quast.report_txt
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
        File? quast_report_txt = FinalizeQuastReportTxt.gcs_path
        File? quast_report_in_various_formats = FinalizeQuastReports.gcs_path

        File? quast_summary_on_all = FinalizeQuastSummaryAll.gcs_path

        File? quast_summary_on_primary = FinalizeQuastIndividualSummary.gcs_path[0]
        File? quast_summary_on_H0 = FinalizeQuastIndividualSummary.gcs_path[1]
        File? quast_summary_on_H1 = FinalizeQuastIndividualSummary.gcs_path[2]
    }
}

task UnzipReference {  
    input {  
        File gzipped_fasta  
    }  

    String output_filename = basename(gzipped_fasta, ".gz")  

    command <<<  
        set -e  
        gunzip -c ~{gzipped_fasta} > ~{output_filename}  
    >>>  

    output {  
        File unzipped_fasta = output_filename  
    }  

    runtime {  
        cpu: 2  
        memory: "16 GiB"  
        disks: "local-disk 20 HDD"  
        docker: "ubuntu:20.04"  
    }  
}  

