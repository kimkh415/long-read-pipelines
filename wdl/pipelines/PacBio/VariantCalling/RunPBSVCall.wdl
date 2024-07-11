version 1.0

import "../../../structs/Structs.wdl"

task Call {
    input {
        Array[File] svsigs
        File ref_fasta
        File ref_fasta_fai
        Boolean ccs
        String prefix
        String zones
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        svsigs:            "per-chromosome *.svsig.gz files"
        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        ccs:               "use optimizations for CCS data"
        prefix:            "prefix for output"
    }

    Int disk_size = 2*ceil(size(svsigs, "GiB") + size([ref_fasta, ref_fasta_fai], "GiB"))

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        pbsv call -j $num_core --log-level INFO ~{true='--ccs' false='' ccs} \
            ~{ref_fasta} \
            ~{sep=' ' svsigs} \
            ~{prefix}.pbsv.pre.vcf

        cat ~{prefix}.pbsv.pre.vcf | grep -v -e '##fileDate' > ~{prefix}.pbsv.vcf
    >>>

    output {
        File vcf = "~{prefix}.pbsv.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             96,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        zones: zones
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


workflow RunPBSVCall {

    meta {
        description: "Run PBSV call"
    }

    parameter_meta {
        svsigs:            "input svsig files output from pbsv discover"
        is_ccs:            "if CCS reads"
        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        prefix:            "prefix for output"
        zones:             "zones to run in"
    }

    input {
        Array[File] svsigs
        Boolean is_ccs

        File ref_fasta
        File ref_fasta_fai
        String prefix

        String zones
    }

    call Call {
        input:
            svsigs        = [ Discover.svsig ],
            ref_fasta     = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ccs           = is_ccs,
            prefix        = prefix,
            zones         = zones
    }

    output {
        File vcf = pbsv.Call.vcf
    }
}


