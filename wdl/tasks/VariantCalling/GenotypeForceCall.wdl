version 1.0


import "../../structs/Structs.wdl"


workflow Sniffles2Genotype {

    meta {
        description: "Use Sniffles2 to call sample genotype on known SVs"
    }

    parameter_meta {
        # input
        inputVCF:        "GCS path to a VCF file with variants to call genotypes on"
        sampleBAM:      "GCS path to an aligned BAM file from a samples"
        sampleBAI:      "GCS path to an aligned BAM file indices from a samples"
        sampleID:       "matching sample ID of the BAM"
        prefix:          "prefix for output file"
        # output
        genotyped_vcf:   "[OUTPUT] vcf output with genotypes"
    }

    input {
        File       inputVCF
        File       sampleBAM
        File       sampleBAI
        String     sampleID
        String     prefix
    }

    call AddGenotype {
        input:
            vcf = inputVCF,
            bam = sampleBAM,
            bai = sampleBAI,
            sample_id = sampleID,
            prefix = prefix
    }

    output {
         File? genotyped_vcf = AddGenotype.geno_vcf
    }
}


task AddGenotype {

    meta {
        description: "This task finds sample genotypes on known SVs based on sample alignment."
    }

    input {
        File vcf
        File bam
        File bai
        String sample_id
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        vcf:              "input VCF of known SVs"
        bam:              "input BAM from which to call genotypes"
        bai:              "index accompanying the BAM"
        sample_id:        "Sample ID"
        prefix:           "prefix"
    }

    Int cpus = 8
    Int disk_size = 2*ceil(size([bam, bai], "GB"))
    String vcf_output = "~{prefix}.~{sample_id}.genotyped.vcf"

    command <<<
        set -eux

        # Install system dependencies
        apt-get update -qq
        apt-get install -y python3-pip build-essential

        # Install latest Sniffles2
        pip install --upgrade pip
        pip install --upgrade sniffles

        # Verify installation and version
        sniffles --version

        sniffles -t ~{cpus} \
                 --sample-id ~{sample_id} \
                 --input ~{bam} \
                 --genotype-vcf ~{vcf} \
                 --vcf ~{vcf_output}
    >>>

    output {
        File geno_vcf = "~{vcf_output}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             46,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "python:3.10-slim"
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
    }
}

