version 1.0

import "../../structs/Structs.wdl"

workflow FragmentReadsRAFT {

    meta {
        description: "After fragmenting the reads using RAFT, we run two HiFiasm jobs, one for getting alternative contigs and one for getting the haplotigs."
    }
    parameter_meta {
        reads:    "reads (in fasta or fastq format, compressed or uncompressed)"
        genome_length:  "estimated length of genome in bases"
        raft_disk_size:  "disk space in GB for RunRAFT"
        est_coverage:   "integer of estimated coverage from raw reads"
        error_correct_read_fa:    "from hifiasm error correct reads"
        overlaps_paf:     "from hifiasm all to all alignment overlaps"
    }

    input {
        File reads
        String genome_length
        Int? est_coverage
        File? error_correct_read_fa
        File? overlaps_paf
        Int raft_disk_size = 300
        String zones = "us-central1-a us-central1-b us-central1-c"
    }

    if (!defined(est_coverage)) {
        call EstimateCoverage {
            input:
                input_fastq = reads,
                estimated_genome_size = genome_length
        }
    }

    if (!defined(error_correct_read_fa)) {
        call GetErrorCorrectedReads {
            input:
                reads = reads,
                zones = zones
        }
    }

    if (!defined(overlaps_paf)) {
        call GetOverlaps {
            input:
                ec_reads = select_first([error_correct_read_fa, GetErrorCorrectedReads.ec_reads]),
                zones = zones
        }
    }

    File ec_reads_output = select_first([error_correct_read_fa, GetErrorCorrectedReads.ec_reads])
    File overlaps_output = select_first([overlaps_paf, GetOverlaps.overlaps])
    Int cov_output = select_first([est_coverage, EstimateCoverage.coverage])

    call InstallAndRunRAFT {
        input:
            error_corrected_reads = ec_reads_output,
            overlaps = overlaps_output,
            coverage = cov_output,
            disk_size = raft_disk_size
    }

    output {
        File ec_reads  = ec_reads_output
        File overlaps = overlaps_output
        Int coverage = cov_output
        File fragmented_reads = RunRAFT.fragmented_reads  
        File raft_log = RunRAFT.log
    }
}

task GetErrorCorrectedReads {
    input {
        File reads
        String zones

        RuntimeAttr? runtime_attr_override
    }

    Int proposed_memory = 4 * ceil(size(reads, "GB"))
    Int memory = if proposed_memory < 96 then 96 else proposed_memory  # this 96 magic number is purely empirical
    Int n = memory / 4  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int disk_size = 10 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        time hifiasm \
            -o errorcorrect \
            -t~{num_cpus} \
            --write-ec \
            ~{reads} \
            2> hifiasm_errorcorrect.log

        tree -h .

    >>>

    output {
        File ec_reads = "errorcorrect.ec.fa"
        File log = "hifiasm_errorcorrect.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/biocontainers/hifiasm:0.20.0--h43eeafb_0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones: zones
    }
}

task GetOverlaps {
    input {
        File ec_reads
        String zones

        RuntimeAttr? runtime_attr_override
    }

    Int proposed_memory = 6 * ceil(size(ec_reads, "GB"))
    Int memory = if proposed_memory < 128 then 128 else proposed_memory  # this 96 magic number is purely empirical
    Int n = memory / 4  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 48 then 48 else num_cpus_proposal

    Int disk_size = 10 * ceil(size(ec_reads, "GB"))

    command <<<
        set -euxo pipefail

        time hifiasm \
            -o getOverlaps \
            -t~{num_cpus} \
            --dbg-ovec \
            ~{ec_reads} \
            2> hifiasm_overlaps.log

        cat getOverlaps.0.ovlp.paf getOverlaps.1.ovlp.paf > overlaps.paf
        
        tree -h .
    >>>

    output {
        File overlaps = "overlaps.paf"
        File log = "hifiasm_overlaps.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/biocontainers/hifiasm:0.20.0--h43eeafb_0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones: zones
    }
}

task EstimateCoverage {
    input {
        File input_fastq
        String estimated_genome_size
        Int cpu = 1
        Int memory_gb = 4
        Int max_retries = 3
    }

    command <<<
        set -euxo pipefail
        total_bases=$(zcat ~{input_fastq} | awk '{if(NR%4==2) sum+=length($0)} END{print sum}')
        echo $((total_bases / ~{estimated_genome_size})) > coverage.txt
    >>>

    output {
        Int coverage = read_int("coverage.txt")
    }

    runtime {
        docker: "ubuntu:latest"
        cpu: cpu
        disks: "local-disk 100 HDD"
        memory: "~{memory_gb} GB"
        maxRetries: max_retries
    }
}

task InstallRAFT {
    input {
        Int cpu = 2
        Int memory_gb = 4
        Int max_retries = 1
    }

    command <<<
        set -euxo pipefail

        apt-get update && apt-get install -y git build-essential
        apt-get install -y libz-dev

        mkdir -p raft_install
        cd raft_install

        git clone https://github.com/at-cg/RAFT.git
        cd RAFT
        make

        mkdir -p $PWD/bin
        mv raft $PWD/bin/

        # Output the full path to the executable
        echo $PWD/bin/raft | tee raft_executable_path.txt
        echo $PWD/bin | tee raft_bin_path.txt
    >>>

    output {
        File raft_executable = "raft_install/RAFT/raft_executable_path.txt"
        String raft_bin_path = read_string("raft_install/RAFT/raft_bin_path.txt")
    }

    runtime {
        docker: "ubuntu:latest"
        cpu: cpu
        disks: "local-disk 20 HDD"
        memory: "~{memory_gb} GB"
        maxRetries: max_retries
        continueOnReturnCode: [0, 1]  # Allow both 0 and 1 as valid return codes
    }
}

task RunRAFT {
    input {
        File error_corrected_reads
        File overlaps
        Int coverage
        String raft_bin_path
        Int cpu = 32
        Int memory_gb = 100
        Int disk_size = 100
        Int max_retries = 3
    }

    command <<<
        set -euxo pipefail
        # Add RAFT to the PATH
        export PATH=$PATH:~{raft_bin_path}

        # Now we can use RAFT directly
        raft -e ~{coverage} -o fragmented ~{error_corrected_reads} ~{overlaps} > raft.log
    >>>

    output {
        File fragmented_reads = "fragmented.reads.fasta"
        File log = "raft.log"
    }

    runtime {
        docker: "ubuntu:latest"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        memory: "~{memory_gb} GB"
        maxRetries: max_retries
    }
}

task InstallAndRunRAFT {
    input {
        File error_corrected_reads
        File overlaps
        Int coverage
        Int cpu = 32  # Using larger value from RunRAFT
        Int memory_gb = 100  # Using larger value from RunRAFT
        Int disk_size = 100
        Int max_retries = 3  # Using larger value from RunRAFT
    }

    command <<<
        set -euxo pipefail

        apt-get update && apt-get install -y git build-essential
        apt-get install -y libz-dev

        mkdir -p raft_install
        cd raft_install

        git clone https://github.com/at-cg/RAFT.git
        cd RAFT
        make

        mkdir -p $PWD/bin
        mv raft $PWD/bin/

        RAFT_PATH=$PWD/bin
        cd ../..  # Return to original directory

        # Add RAFT to path
        export PATH=$PATH:$RAFT_PATH

        # Run RAFT
        raft -e ~{coverage} -o fragmented ~{error_corrected_reads} ~{overlaps} > raft.log
    >>>  

    output {  
        File fragmented_reads = "fragmented.reads.fasta"  
        File log = "raft.log"  
    }  

    runtime {  
        docker: "ubuntu:latest"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        memory: "~{memory_gb} GB"
        maxRetries: max_retries
    }
}
