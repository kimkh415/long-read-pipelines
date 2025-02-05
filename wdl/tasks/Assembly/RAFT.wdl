version 1.0

task InstallRAFT {
    input {
        Int cpu = 2
        Int memory_gb = 4
        Int max_retries = 3
    }

    command <<<
        set -e
        mkdir -p raft_install
        cd raft_install

        git clone https://github.com/at-cg/RAFT.git
        cd RAFT
        make

        mkdir -p $PWD/bin
        mv raft $PWD/bin/

        # Output the full path to the executable
        echo $PWD/bin/raft > raft_executable_path.txt
        echo $PWD/bin > raft_bin_path.txt

        # Verify installation
        $PWD/bin/raft --help
    >>>

    output {
        File raft_executable = read_string("raft_executable_path.txt")
        String raft_bin_path = read_string("raft_bin_path.txt")
    }

    runtime {
        docker: "ubuntu:latest"
        cpu: cpu
        memory: "~{memory_gb} GB"
        maxRetries: max_retries
    }
}

task EstimateCoverage {
    input {
        File input_fastq
        Int estimated_genome_size
        Int cpu = 1
        Int memory_gb = 4
        Int max_retries = 3
    }

    command <<<
        set -e
        total_bases=$(zcat ~{input_fastq} | awk '{if(NR%4==2) sum+=length($0)} END{print sum}')
        echo $((total_bases / ~{estimated_genome_size})) > coverage.txt
    >>>

    output {
        Int coverage = read_int("coverage.txt")
    }

    runtime {
        docker: "ubuntu:latest"
        cpu: cpu
        memory: "~{memory_gb} GB"
        maxRetries: max_retries
    }
}

task MergeOverlaps {
    input {
        File cis_overlaps
        File trans_overlaps
        Int cpu = 1
        Int memory_gb = 4
        Int max_retries = 3
    }

    command <<<
        set -e
        cat ~{cis_overlaps} ~{trans_overlaps} > merged_overlaps.paf
    >>>

    output {
        File merged_overlaps = "merged_overlaps.paf"
    }

    runtime {
        docker: "ubuntu:latest"
        cpu: cpu
        memory: "~{memory_gb} GB"
        maxRetries: max_retries
    }
}

task RunRAFT {
    input {
        File error_corrected_reads
        File overlaps
        Int coverage
        String raft_bin_path
        Int cpu = 4
        Int memory_gb = 16
        Int max_retries = 3
    }

    command <<<
        set -e
        # Add RAFT to the PATH
        export PATH=$PATH:~{raft_bin_path}

        # Now we can use RAFT directly
        raft -e ~{coverage} -o fragmented ~{error_corrected_reads} ~{overlaps}
    >>>

    output {
        File fragmented_reads = "fragmented.reads.fasta"
    }

    runtime {
        docker: "ubuntu:latest"
        cpu: cpu
        memory: "~{memory_gb} GB"
        maxRetries: max_retries
    }
}
