version 1.0

struct RuntimeAttributes {
    Int? cpu
    Int? command_mem_gb
    Int? additional_mem_gb
    Int? disk_size_gb
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible
    Int? max_retries
}

workflow KAGE2Index {
    input {
        File panel_vcf_gz
        File panel_vcf_gz_tbi
        File reference_fasta
        File reference_fasta_fai
        Array[String] chromosomes
        String output_prefix

        String kage2_docker
        File? monitoring_script
        RuntimeAttributes? kage2_runtime_attributes
    }

    call Index {
        input:
            panel_vcf_gz = panel_vcf_gz,
            panel_vcf_gz_tbi = panel_vcf_gz_tbi,
            reference_fasta = reference_fasta,
            chromosomes = chromosomes,
            output_prefix = output_prefix,
            docker = kage2_docker,
            monitoring_script = monitoring_script,
            runtime_attributes = kage2_runtime_attributes
    }

    output {
        File kage2_index = Index.kage2_index
    }
}

task Index {
    input {
        File panel_vcf_gz
        File panel_vcf_gz_tbi
        File reference_fasta
        Array[String] chromosomes
        String output_prefix

        String docker
        File? monitoring_script
        Int? kmer_length = 31
        String? extra_args = "--no-helper-model True"

        RuntimeAttributes runtime_attributes = {}
    }

    Int num_chromosomes = length(chromosomes)

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # subset reference to chromosomes and replace non-ACGTN with N
        samtools faidx -r <(echo -e "~{sep="\n" chromosomes}") ~{reference_fasta} | \
            sed -e '/>/!s/[^ACTGN]/N/g' > reference.subset.fa

        # subset panel VCF to chromosomes, split to biallelic, fill AF
        # TODO do we need to sort?
        bcftools norm -m -any ~{panel_vcf_gz} -r ~{sep="," chromosomes} | \
            bcftools +fill-tags -O panel.subset.vcf -- -t AF

        NPROC=$(nproc)
        NUM_THREADS=$NPROC

        kage index \
            -o ~{output_prefix}.index \
            -r reference.subset.fa \
            -t $NUM_THREADS \
            -v panel.subset.vcf \
            -k ~{kmer_length} \
            ~{extra_args}
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File kage2_index = "~{output_prefix}.index"
    }
}
